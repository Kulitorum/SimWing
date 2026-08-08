#include "canonical_case.h"
#include "open_piston_case.h"
#include "piston_case.h"
#include "viewer_protocol.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <spawn.h>
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

extern char** environ;
#endif

namespace {

constexpr std::uint64_t defaultSteps = 600;
constexpr std::uint64_t maximumSteps = 10'000'000;

enum class WorkerCase {
    Structural,
    Piston,
    OpenPiston,
};

struct Options {
    std::uint64_t steps = defaultSteps;
    std::filesystem::path tracePath;
    bool viewer = true;
    bool help = false;
    WorkerCase workerCase = WorkerCase::Structural;
};

void printUsage(FILE* stream) {
    std::fprintf(
        stream,
        "Usage: simwing-fsi [--case structural|piston|open-piston] [--steps N]\n"
        "                   [--trace PATH]\n"
        "                   [--viewer|--no-viewer]\n"
        "\n"
        "Runs a canonical Qt-free numerical case and writes a completed diagnostic\n"
        "trace. 'structural' is the original analytic XPBD harness; 'piston' runs\n"
        "the face-resolved fluid -> transfer -> temporal coupling -> XPBD path;\n"
        "'open-piston' adds connected-fluid pressure reaction, partial-cell motion,\n"
        "an independently closed opening-flux GCL ledger, and exact one-face\n"
        "topology rebasing; consuming its resolved opening is rejected. Interactive\n"
        "runs launch the sibling simwing-viewer with --follow; --no-viewer is\n"
        "unthrottled for tests and CI.\n");
}

bool parseWorkerCase(const std::string_view text, WorkerCase& workerCase) {
    if (text == "structural") {
        workerCase = WorkerCase::Structural;
        return true;
    }
    if (text == "piston") {
        workerCase = WorkerCase::Piston;
        return true;
    }
    if (text == "open-piston") {
        workerCase = WorkerCase::OpenPiston;
        return true;
    }
    return false;
}

bool parseUnsigned(std::string_view text, std::uint64_t& value) {
    if (text.empty()) {
        return false;
    }
    const char* begin = text.data();
    const char* end = begin + text.size();
    const auto [position, error] = std::from_chars(begin, end, value);
    return error == std::errc{} && position == end;
}

bool parseOptions(int argc,
                  char* argv[],
                  Options& options,
                  std::string& error) {
    bool viewerRequested = false;
    bool noViewerRequested = false;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            options.help = true;
        } else if (argument == "--viewer") {
            viewerRequested = true;
            options.viewer = true;
        } else if (argument == "--no-viewer") {
            noViewerRequested = true;
            options.viewer = false;
        } else if (argument == "--case") {
            if (++index >= argc
                || !parseWorkerCase(argv[index], options.workerCase)) {
                error = "--case requires 'structural', 'piston', or 'open-piston'";
                return false;
            }
        } else if (argument.starts_with("--case=")) {
            if (!parseWorkerCase(argument.substr(7), options.workerCase)) {
                error = "--case requires 'structural', 'piston', or 'open-piston'";
                return false;
            }
        } else if (argument == "--steps") {
            if (++index >= argc
                || !parseUnsigned(argv[index], options.steps)) {
                error = "--steps requires an unsigned integer";
                return false;
            }
        } else if (argument.starts_with("--steps=")) {
            if (!parseUnsigned(argument.substr(8), options.steps)) {
                error = "--steps requires an unsigned integer";
                return false;
            }
        } else if (argument == "--trace") {
            if (++index >= argc || argv[index][0] == '\0') {
                error = "--trace requires a path";
                return false;
            }
            options.tracePath = std::filesystem::path(argv[index]);
        } else if (argument.starts_with("--trace=")) {
            if (argument.size() == 8) {
                error = "--trace requires a path";
                return false;
            }
            options.tracePath = std::filesystem::path(argument.substr(8));
        } else {
            error = "unknown option: " + std::string(argument);
            return false;
        }
    }
    if (viewerRequested && noViewerRequested) {
        error = "--viewer and --no-viewer are mutually exclusive";
        return false;
    }
    if (options.steps > maximumSteps) {
        error = "--steps exceeds the worker safety limit";
        return false;
    }
    return true;
}

std::uint64_t processId() {
#ifdef _WIN32
    return static_cast<std::uint64_t>(GetCurrentProcessId());
#else
    return static_cast<std::uint64_t>(getpid());
#endif
}

std::filesystem::path defaultTracePath() {
    return std::filesystem::temp_directory_path()
        / ("simwing-fsi-" + std::to_string(processId()) + ".swtrace");
}

std::filesystem::path siblingViewerPath(const char* workerArgumentZero) {
    std::error_code error;
#ifdef _WIN32
    std::filesystem::path worker;
    for (std::size_t size = 512; size <= 32768; size *= 2) {
        std::vector<wchar_t> buffer(size);
        const DWORD copied = GetModuleFileNameW(
            nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
        if (copied > 0 && copied < buffer.size()) {
            worker = std::filesystem::path(
                std::wstring_view(buffer.data(), copied));
            break;
        }
    }
    if (worker.empty()) {
        worker = std::filesystem::absolute(
            std::filesystem::path(workerArgumentZero), error);
    }
#else
    std::filesystem::path worker;
#if defined(__linux__)
    worker = std::filesystem::read_symlink("/proc/self/exe", error);
#elif defined(__APPLE__)
    std::uint32_t size = 0;
    static_cast<void>(_NSGetExecutablePath(nullptr, &size));
    std::vector<char> buffer(size);
    if (size > 0 && _NSGetExecutablePath(buffer.data(), &size) == 0) {
        worker = std::filesystem::path(buffer.data());
    }
#endif
    if (worker.empty() || error) {
        error.clear();
        worker = std::filesystem::absolute(
            std::filesystem::path(workerArgumentZero), error);
    }
#endif
    if (error) {
        worker = std::filesystem::path(workerArgumentZero);
    }
    worker = std::filesystem::weakly_canonical(worker, error);
#ifdef _WIN32
    constexpr const char* viewerName = "simwing-viewer.exe";
#else
    constexpr const char* viewerName = "simwing-viewer";
#endif
    return worker.parent_path() / viewerName;
}

#ifdef _WIN32
std::wstring quoteWindowsArgument(const std::wstring& argument) {
    std::wstring quoted;
    quoted.push_back(L'"');
    std::size_t backslashes = 0;
    for (const wchar_t character : argument) {
        if (character == L'\\') {
            ++backslashes;
        } else if (character == L'"') {
            quoted.append(backslashes * 2 + 1, L'\\');
            quoted.push_back(L'"');
            backslashes = 0;
        } else {
            quoted.append(backslashes, L'\\');
            backslashes = 0;
            quoted.push_back(character);
        }
    }
    quoted.append(backslashes * 2, L'\\');
    quoted.push_back(L'"');
    return quoted;
}
#endif

bool launchViewer(const std::filesystem::path& viewer,
                  const std::filesystem::path& trace,
                  std::string& error) {
    if (!std::filesystem::is_regular_file(viewer)) {
        error = "viewer executable not found beside worker: "
            + viewer.string();
        return false;
    }
#ifdef _WIN32
    std::wstring commandLine = quoteWindowsArgument(viewer.native())
        + L" --follow " + quoteWindowsArgument(trace.native());
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process{};
    if (!CreateProcessW(viewer.c_str(), commandLine.data(), nullptr, nullptr,
                        FALSE, CREATE_NEW_PROCESS_GROUP, nullptr,
                        viewer.parent_path().c_str(), &startup, &process)) {
        error = "CreateProcessW failed with error "
            + std::to_string(GetLastError());
        return false;
    }
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    return true;
#else
    std::string executable = viewer.string();
    std::string followArgument = "--follow";
    std::string traceArgument = trace.string();
    std::vector<char*> arguments{
        executable.data(),
        followArgument.data(),
        traceArgument.data(),
        nullptr,
    };
    pid_t child = 0;
    const int result = posix_spawn(&child, executable.c_str(), nullptr,
                                   nullptr, arguments.data(), environ);
    if (result != 0) {
        error = "posix_spawn failed with error " + std::to_string(result);
        return false;
    }
    return true;
#endif
}

bool flushTrace(std::ofstream& output, std::string& error) {
    output.flush();
    if (!output) {
        error = "failed to flush completed trace record";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char* argv[]) {
    Options options;
    std::string error;
    if (!parseOptions(argc, argv, options, error)) {
        std::fprintf(stderr, "%s\n", error.c_str());
        printUsage(stderr);
        return 2;
    }
    if (options.help) {
        printUsage(stdout);
        return 0;
    }

    try {
        if (options.tracePath.empty()) {
            options.tracePath = defaultTracePath();
        }
        options.tracePath = std::filesystem::absolute(options.tracePath);

        std::ofstream output(options.tracePath,
                             std::ios::binary | std::ios::trunc);
        if (!output) {
            std::fprintf(stderr, "cannot open trace for writing: %s\n",
                         options.tracePath.string().c_str());
            return 1;
        }

        const auto run = [&](auto& simulation) -> int {
            simwing::viewer::TraceWriter writer(output);
            if (!writer.writeHeader(simulation.traceHeader())
                || !flushTrace(output, error)) {
                if (error.empty()) {
                    error = writer.error().message;
                }
                std::fprintf(stderr, "trace header failed: %s\n",
                             error.c_str());
                return 1;
            }

            if (options.viewer) {
                const std::filesystem::path viewer =
                    siblingViewerPath(argv[0]);
                if (!launchViewer(viewer, options.tracePath, error)) {
                    std::fprintf(stderr, "%s\n", error.c_str());
                    return 1;
                }
            }

            auto nextFrameTime = std::chrono::steady_clock::now();
            for (std::uint64_t step = 0; step < options.steps; ++step) {
                const simwing::viewer::DiagnosticFrame frame =
                    simulation.advance();
                if (!writer.writeFrame(frame)
                    || !flushTrace(output, error)) {
                    if (error.empty()) {
                        error = writer.error().message;
                    }
                    std::fprintf(stderr, "trace frame failed: %s\n",
                                 error.c_str());
                    return 1;
                }
                if (options.viewer) {
                    nextFrameTime += std::chrono::duration_cast<
                        std::chrono::steady_clock::duration>(
                        std::chrono::duration<double>(
                            simulation.stepSettings().timeStepSeconds));
                    std::this_thread::sleep_until(nextFrameTime);
                }
            }

            if (!writer.finish() || !flushTrace(output, error)) {
                if (error.empty()) {
                    error = writer.error().message;
                }
                std::fprintf(stderr, "trace completion failed: %s\n",
                             error.c_str());
                return 1;
            }

            const simwing::fsi::StructureCheckpoint checkpoint =
                simulation.structure().checkpoint();
            simwing::fsi::StructureVector3 minimum =
                checkpoint.nodes.front().positionMeters;
            simwing::fsi::StructureVector3 maximum = minimum;
            for (const simwing::fsi::StructureNodeState& node
                 : checkpoint.nodes) {
                minimum.x = std::min(minimum.x, node.positionMeters.x);
                minimum.y = std::min(minimum.y, node.positionMeters.y);
                minimum.z = std::min(minimum.z, node.positionMeters.z);
                maximum.x = std::max(maximum.x, node.positionMeters.x);
                maximum.y = std::max(maximum.y, node.positionMeters.y);
                maximum.z = std::max(maximum.z, node.positionMeters.z);
            }
            const simwing::fsi::StructureDiagnostics diagnostics =
                simulation.structure().diagnostics();
            std::printf("simwing-fsi completed %llu step(s), t=%.9g s, "
                        "bounds=[%.6g %.6g %.6g]-[%.6g %.6g %.6g] m, "
                        "max-strain=%.6g, trace=%s\n",
                        static_cast<unsigned long long>(
                            checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        minimum.x, minimum.y, minimum.z,
                        maximum.x, maximum.y, maximum.z,
                        diagnostics.maximumAbsoluteMembraneStrain,
                        options.tracePath.string().c_str());
            return 0;
        };

        if (options.workerCase == WorkerCase::Piston) {
            simwing::fsi::CoupledPistonCase simulation;
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::OpenPiston) {
            simwing::fsi::OpenPistonCase simulation;
            return run(simulation);
        }
        simwing::fsi::CanonicalStructuralCase simulation;
        return run(simulation);
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "simwing-fsi failed: %s\n", exception.what());
        return 1;
    }
}
