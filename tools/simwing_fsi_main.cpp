#include "canonical_case.h"
#include "open_piston_case.h"
#include "periodic_flow_case.h"
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
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <type_traits>
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
    PeriodicFlow,
};

struct Options {
    std::uint64_t steps = defaultSteps;
    std::filesystem::path tracePath;
    std::filesystem::path checkpointInputPath;
    std::filesystem::path checkpointOutputPath;
    std::uint64_t checkpointEvery = 0;
    bool viewer = true;
    bool help = false;
    WorkerCase workerCase = WorkerCase::Structural;
};

void printUsage(FILE* stream) {
    std::fprintf(
        stream,
        "Usage: simwing-fsi [--case structural|piston|open-piston|periodic-flow]\n"
        "                   [--steps N]\n"
        "                   [--trace PATH]\n"
        "                   [--checkpoint-in PATH]\n"
        "                   [--checkpoint-out PATH]\n"
        "                   [--checkpoint-every N]\n"
        "                   [--viewer|--no-viewer]\n"
        "\n"
        "Runs a canonical Qt-free numerical case and writes a completed diagnostic\n"
        "trace. 'structural' is the original analytic XPBD harness; 'piston' runs\n"
        "the face-resolved fluid -> transfer -> temporal coupling -> XPBD path;\n"
        "'open-piston' adds connected-fluid pressure reaction, partial-cell motion,\n"
        "an independently closed opening-flux GCL ledger, and exact one-face\n"
        "topology rebasing; consuming its resolved opening is rejected.\n"
        "'periodic-flow' advances the bounded Strang/SSPRK2 Taylor-Green CFD\n"
        "canonical and publishes cell-centred pressure/velocity points. Its optional\n"
        "checkpoint paths restore/save exact accepted state; --checkpoint-every\n"
        "autosaves at absolute accepted-step multiples and the final state. --steps\n"
        "counts additional intervals. Interactive runs launch the sibling viewer\n"
        "with --follow; --no-viewer is unthrottled for tests and CI.\n");
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
    if (text == "periodic-flow") {
        workerCase = WorkerCase::PeriodicFlow;
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
                error = "--case requires 'structural', 'piston', "
                    "'open-piston', or 'periodic-flow'";
                return false;
            }
        } else if (argument.starts_with("--case=")) {
            if (!parseWorkerCase(argument.substr(7), options.workerCase)) {
                error = "--case requires 'structural', 'piston', "
                    "'open-piston', or 'periodic-flow'";
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
        } else if (argument == "--checkpoint-in") {
            if (++index >= argc || argv[index][0] == '\0') {
                error = "--checkpoint-in requires a path";
                return false;
            }
            options.checkpointInputPath = std::filesystem::path(argv[index]);
        } else if (argument.starts_with("--checkpoint-in=")) {
            if (argument.size() == 16) {
                error = "--checkpoint-in requires a path";
                return false;
            }
            options.checkpointInputPath =
                std::filesystem::path(argument.substr(16));
        } else if (argument == "--checkpoint-out") {
            if (++index >= argc || argv[index][0] == '\0') {
                error = "--checkpoint-out requires a path";
                return false;
            }
            options.checkpointOutputPath =
                std::filesystem::path(argv[index]);
        } else if (argument.starts_with("--checkpoint-out=")) {
            if (argument.size() == 17) {
                error = "--checkpoint-out requires a path";
                return false;
            }
            options.checkpointOutputPath =
                std::filesystem::path(argument.substr(17));
        } else if (argument == "--checkpoint-every") {
            if (++index >= argc
                || !parseUnsigned(argv[index], options.checkpointEvery)
                || options.checkpointEvery == 0) {
                error = "--checkpoint-every requires a positive integer";
                return false;
            }
        } else if (argument.starts_with("--checkpoint-every=")) {
            if (!parseUnsigned(
                    argument.substr(19), options.checkpointEvery)
                || options.checkpointEvery == 0) {
                error = "--checkpoint-every requires a positive integer";
                return false;
            }
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
    if ((!options.checkpointInputPath.empty()
         || !options.checkpointOutputPath.empty())
        && options.workerCase != WorkerCase::PeriodicFlow) {
        error = "checkpoint paths require --case periodic-flow";
        return false;
    }
    if (options.checkpointEvery != 0
        && options.checkpointOutputPath.empty()) {
        error = "--checkpoint-every requires --checkpoint-out";
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

bool samePath(const std::filesystem::path& first,
              const std::filesystem::path& second) {
    if (first == second) {
        return true;
    }
#ifdef _WIN32
    if (CompareStringOrdinal(
            first.c_str(), -1, second.c_str(), -1, TRUE) == CSTR_EQUAL) {
        return true;
    }
#endif
    std::error_code error;
    const bool equivalent = std::filesystem::exists(first, error)
        && !error
        && std::filesystem::exists(second, error)
        && !error
        && std::filesystem::equivalent(first, second, error);
    return equivalent && !error;
}

std::filesystem::path normalizedAbsolutePath(
    const std::filesystem::path& path) {
    std::filesystem::path result =
        std::filesystem::absolute(path).lexically_normal();
    std::error_code error;
    const std::filesystem::path canonical =
        std::filesystem::weakly_canonical(result, error);
    return error ? result : canonical;
}

std::filesystem::path checkpointTemporaryPath(
    const std::filesystem::path& path) {
    std::filesystem::path temporary = path;
    temporary += ".tmp-" + std::to_string(processId());
    return temporary;
}

bool readPeriodicFlowCheckpoint(
    const std::filesystem::path& path,
    simwing::fsi::PeriodicFlowCaseCheckpoint& checkpoint,
    std::string& error) {
    const simwing::fsi::PeriodicFlowCaseCheckpointLimits limits;
    std::error_code fileError;
    const std::uintmax_t size = std::filesystem::file_size(path, fileError);
    if (fileError) {
        error = "cannot inspect periodic-flow checkpoint: " + path.string();
        return false;
    }
    if (size > limits.maximumBytes
        || size > std::numeric_limits<std::size_t>::max()) {
        error = "periodic-flow checkpoint exceeds the byte limit";
        return false;
    }
    std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        error = "cannot open periodic-flow checkpoint: " + path.string();
        return false;
    }
    if (!bytes.empty()) {
        input.read(
            reinterpret_cast<char*>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
    }
    if (!input || input.peek() != std::char_traits<char>::eof()) {
        error = "failed to read complete periodic-flow checkpoint";
        return false;
    }
    simwing::fsi::PeriodicFlowCaseCheckpointError protocolError;
    if (!simwing::fsi::deserializePeriodicFlowCaseCheckpoint(
            bytes, checkpoint, &protocolError, limits)) {
        error = protocolError.message;
        return false;
    }
    return true;
}

bool writePeriodicFlowCheckpoint(
    const std::filesystem::path& path,
    const simwing::fsi::PeriodicFlowCaseCheckpoint& checkpoint,
    std::string& error) {
    std::vector<std::uint8_t> bytes;
    simwing::fsi::PeriodicFlowCaseCheckpointError protocolError;
    if (!simwing::fsi::serializePeriodicFlowCaseCheckpoint(
            checkpoint, bytes, &protocolError)) {
        error = protocolError.message;
        return false;
    }
    const std::filesystem::path temporary = checkpointTemporaryPath(path);
    std::error_code temporaryError;
    if (std::filesystem::exists(temporary, temporaryError)
        || temporaryError) {
        error = "refusing to overwrite temporary checkpoint path: "
            + temporary.string();
        return false;
    }
    {
        std::ofstream output(
            temporary, std::ios::binary | std::ios::trunc);
        if (!output) {
            error = "cannot open temporary checkpoint for writing: "
                + temporary.string();
            return false;
        }
        output.write(
            reinterpret_cast<const char*>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
        output.flush();
        if (!output) {
            error = "failed to write complete periodic-flow checkpoint";
            output.close();
            std::error_code removeError;
            std::filesystem::remove(temporary, removeError);
            return false;
        }
    }

#ifdef _WIN32
    if (!MoveFileExW(
            temporary.c_str(), path.c_str(),
            MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) {
        error = "cannot atomically replace periodic-flow checkpoint, error "
            + std::to_string(GetLastError());
        std::error_code removeError;
        std::filesystem::remove(temporary, removeError);
        return false;
    }
#else
    std::error_code renameError;
    std::filesystem::rename(temporary, path, renameError);
    if (renameError) {
        error = "cannot atomically replace periodic-flow checkpoint: "
            + renameError.message();
        std::error_code removeError;
        std::filesystem::remove(temporary, removeError);
        return false;
    }
#endif
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
        options.tracePath = normalizedAbsolutePath(options.tracePath);
        if (!options.checkpointInputPath.empty()) {
            options.checkpointInputPath = normalizedAbsolutePath(
                options.checkpointInputPath);
        }
        if (!options.checkpointOutputPath.empty()) {
            options.checkpointOutputPath = normalizedAbsolutePath(
                options.checkpointOutputPath);
        }
        if ((!options.checkpointInputPath.empty()
             && samePath(options.tracePath, options.checkpointInputPath))
            || (!options.checkpointOutputPath.empty()
                && samePath(
                    options.tracePath, options.checkpointOutputPath))
            || (!options.checkpointOutputPath.empty()
                && samePath(
                    options.tracePath,
                    checkpointTemporaryPath(
                        options.checkpointOutputPath)))) {
            std::fprintf(
                stderr,
                "trace and checkpoint paths must be different files\n");
            return 2;
        }

        std::optional<simwing::fsi::PeriodicFlowCaseCheckpoint>
            restoredPeriodicFlowCheckpoint;
        if (!options.checkpointInputPath.empty()) {
            restoredPeriodicFlowCheckpoint.emplace();
            if (!readPeriodicFlowCheckpoint(
                    options.checkpointInputPath,
                    *restoredPeriodicFlowCheckpoint,
                    error)) {
                std::fprintf(stderr, "checkpoint restore failed: %s\n",
                             error.c_str());
                return 1;
            }
        }

        std::ofstream output(options.tracePath,
                             std::ios::binary | std::ios::trunc);
        if (!output) {
            std::fprintf(stderr, "cannot open trace for writing: %s\n",
                         options.tracePath.string().c_str());
            return 1;
        }

        const auto run = [&](auto& simulation) -> int {
            std::uint64_t lastCheckpointStep =
                std::numeric_limits<std::uint64_t>::max();
            std::uint64_t checkpointWriteCount = 0;
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

            const double stepSeconds = [&] {
                if constexpr (requires {
                                  simulation.stepSettings()
                                      .flow.timeStepSeconds;
                              }) {
                    return simulation.stepSettings().flow.timeStepSeconds;
                } else {
                    return simulation.stepSettings().timeStepSeconds;
                }
            }();
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
                if constexpr (std::is_same_v<
                                  std::remove_cvref_t<decltype(simulation)>,
                                  simwing::fsi::PeriodicFlowCase>) {
                    if (options.checkpointEvery != 0
                        && simulation.acceptedStepCount()
                               % options.checkpointEvery == 0) {
                        if (!writePeriodicFlowCheckpoint(
                                options.checkpointOutputPath,
                                simulation.checkpoint(),
                                error)) {
                            std::fprintf(
                                stderr,
                                "checkpoint autosave failed: %s\n",
                                error.c_str());
                            return 1;
                        }
                        lastCheckpointStep =
                            simulation.acceptedStepCount();
                        ++checkpointWriteCount;
                    }
                }
                if (options.viewer) {
                    nextFrameTime += std::chrono::duration_cast<
                        std::chrono::steady_clock::duration>(
                        std::chrono::duration<double>(stepSeconds));
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

            if constexpr (std::is_same_v<
                              std::remove_cvref_t<decltype(simulation)>,
                              simwing::fsi::PeriodicFlowCase>) {
                if (!options.checkpointOutputPath.empty()
                    && lastCheckpointStep
                        != simulation.acceptedStepCount()) {
                    if (!writePeriodicFlowCheckpoint(
                            options.checkpointOutputPath,
                            simulation.checkpoint(),
                            error)) {
                        std::fprintf(stderr, "checkpoint save failed: %s\n",
                                     error.c_str());
                        return 1;
                    }
                    lastCheckpointStep = simulation.acceptedStepCount();
                    ++checkpointWriteCount;
                }
            }

            if constexpr (requires { simulation.structure(); }) {
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
                std::printf(
                    "simwing-fsi completed %llu step(s), t=%.9g s, "
                    "bounds=[%.6g %.6g %.6g]-[%.6g %.6g %.6g] m, "
                    "max-strain=%.6g, trace=%s\n",
                    static_cast<unsigned long long>(
                        checkpoint.acceptedStepCount),
                    checkpoint.simulationTimeSeconds,
                    minimum.x, minimum.y, minimum.z,
                    maximum.x, maximum.y, maximum.z,
                    diagnostics.maximumAbsoluteMembraneStrain,
                    options.tracePath.string().c_str());
            } else {
                const auto& diagnostics = simulation.diagnostics();
                std::printf(
                    "simwing-fsi completed %llu periodic-flow step(s), "
                    "t=%.9g s, substeps=%zu, max-CFL=%.6g, "
                    "divergence-L2=%.6g 1/s, energy=%.9g J, "
                    "checkpoint-writes=%llu, trace=%s\n",
                    static_cast<unsigned long long>(
                        simulation.acceptedStepCount()),
                    simulation.simulationTimeSeconds(),
                    diagnostics.plannedSubstepCount,
                    diagnostics.maximumObservedOutgoingCourantNumber,
                    diagnostics.finalDivergenceL2PerSecond,
                    diagnostics.kineticEnergyAfterJoules,
                    static_cast<unsigned long long>(checkpointWriteCount),
                    options.tracePath.string().c_str());
            }
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
        if (options.workerCase == WorkerCase::PeriodicFlow) {
            simwing::fsi::PeriodicFlowCase simulation;
            if (restoredPeriodicFlowCheckpoint) {
                simulation.restore(*restoredPeriodicFlowCheckpoint);
            }
            return run(simulation);
        }
        simwing::fsi::CanonicalStructuralCase simulation;
        return run(simulation);
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "simwing-fsi failed: %s\n", exception.what());
        return 1;
    }
}
