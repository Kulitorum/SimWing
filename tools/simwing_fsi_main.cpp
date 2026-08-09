#include "canonical_case.h"
#include "hemisphere_case.h"
#include "moving_porous_flow_case.h"
#include "moving_porous_flow_checkpoint_persistence.h"
#include "moving_porous_flow_control.h"
#include "open_piston_case.h"
#include "open_piston_checkpoint_persistence.h"
#include "open_piston_control.h"
#include "periodic_flow_control.h"
#include "periodic_flow_case.h"
#include "piston_case.h"
#include "porous_flow_case.h"
#include "porous_sheet_case.h"
#include "porous_sheet_checkpoint_persistence.h"
#include "porous_sheet_control.h"
#include "pressure_jump_case.h"
#include "projected_flag_case.h"
#include "ram_air_cell_case.h"
#include "scene_pressure_cell_case.h"
#include "scene_pressure_cell_checkpoint_persistence.h"
#include "strong_piston_checkpoint_persistence.h"
#include "strong_piston_control.h"
#include "viewer_protocol.h"
#include "worker_control_stream.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
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
    Hemisphere,
    ProjectedFlag,
    RamAirCell,
    ScenePressureCell,
    Piston,
    StrongPiston,
    OpenPiston,
    PeriodicFlow,
    PorousFlow,
    MovingPorousFlow,
    PorousSheet,
    PressureJump,
};

struct Options {
    std::uint64_t steps = defaultSteps;
    std::filesystem::path tracePath;
    std::filesystem::path checkpointInputPath;
    std::filesystem::path checkpointOutputPath;
    std::uint64_t checkpointEvery = 0;
    bool viewer = true;
    bool controlStdio = false;
    bool mimeticPressureAudit = false;
    bool help = false;
    WorkerCase workerCase = WorkerCase::Structural;
};

void printUsage(FILE* stream) {
    std::fprintf(
        stream,
        "Usage: simwing-fsi [--case structural|hemisphere|flag|ram-cell|pressure-cell|piston|strong-piston|open-piston|periodic-flow|porous-flow|moving-porous-flow|porous-sheet|pressure-jump]\n"
        "                   [--steps N]\n"
        "                   [--trace PATH]\n"
        "                   [--checkpoint-in PATH]\n"
        "                   [--checkpoint-out PATH]\n"
        "                   [--checkpoint-every N]\n"
        "                   [--mimetic-pressure-audit]\n"
        "                   [--control-stdio]\n"
        "                   [--viewer|--no-viewer]\n"
        "\n"
        "Runs a canonical Qt-free numerical case and writes a completed diagnostic\n"
        "trace. 'structural' is the original analytic XPBD harness; 'hemisphere'\n"
        "runs a soft three-point fabric dome under an alternating pressure mode;\n"
        "'flag' maps the complete reaction from a fixed-reference projected gust\n"
        "onto a one-edge-clamped XPBD fabric panel;\n"
        "'ram-cell' maps five fixed-reference cavity-wall reactions onto one\n"
        "shared-node open XPBD fabric cell;\n"
        "'pressure-cell' strongly couples scene-v2 moving cut-volume pressure\n"
        "back to an open XPBD cell while a prescribed mean-flow pump and\n"
        "symmetric viscous/projected nonlinear flow advance its bulk MAC predictor;\n"
        "--mimetic-pressure-audit additionally runs the mixed-hybrid pressure\n"
        "path as a read-only pressure-cell shadow after graph convergence; it\n"
        "reports graph-versus-shadow pressure/load deltas and persists its\n"
        "compact accepted state in pressure-cell checkpoints;\n"
        "'piston' runs\n"
        "the face-resolved fluid -> transfer -> temporal coupling -> XPBD path;\n"
        "'strong-piston' strongly iterates that chain for a light added-mass plate;\n"
        "'open-piston' adds connected-fluid pressure reaction, partial-cell motion,\n"
        "an independently closed opening-flux GCL ledger, and exact one-face\n"
        "topology rebasing; consuming its resolved opening is rejected.\n"
        "'periodic-flow' advances the bounded Strang/SSPRK2 Taylor-Green CFD\n"
        "canonical and publishes cell-centred pressure/velocity points.\n"
        "'porous-flow' advances a pressure-driven Darcy-Forchheimer plug and\n"
        "publishes its porous and periodic gauge-closure planes.\n"
        "'moving-porous-flow' advances a prescribed translating porous plane\n"
        "through both symmetric pressure stages and periodic topology wraps.\n"
        "'porous-sheet' couples a translating Darcy sheet to XPBD while a\n"
        "periodic pump drives fluid through it across one MAC-face rebase.\n"
        "'pressure-jump' repeatedly verifies a static split-region slab and\n"
        "publishes each ordered sharp-interface layer. Moving-porous-flow,\n"
        "pressure-cell, strong-piston, porous-sheet, open-piston, and periodic-flow checkpoint\n"
        "paths\n"
        "restore/save exact accepted state;\n"
        "--checkpoint-every\n"
        "autosaves at absolute accepted-step multiples and the final state. --steps\n"
        "counts additional intervals. --control-stdio instead exchanges bounded\n"
        "binary strong-piston, moving-porous-flow, porous-sheet, open-piston, or\n"
        "periodic-flow\n"
        "commands/responses\n"
        "on stdin/stdout;\n"
        "it rejects\n"
        "--steps, --checkpoint-every, and viewer launch. Interactive runs launch\n"
        "the sibling viewer with --follow; --no-viewer is unthrottled for tests\n"
        "and CI.\n");
}

bool parseWorkerCase(const std::string_view text, WorkerCase& workerCase) {
    if (text == "structural") {
        workerCase = WorkerCase::Structural;
        return true;
    }
    if (text == "hemisphere") {
        workerCase = WorkerCase::Hemisphere;
        return true;
    }
    if (text == "flag") {
        workerCase = WorkerCase::ProjectedFlag;
        return true;
    }
    if (text == "ram-cell") {
        workerCase = WorkerCase::RamAirCell;
        return true;
    }
    if (text == "pressure-cell") {
        workerCase = WorkerCase::ScenePressureCell;
        return true;
    }
    if (text == "piston") {
        workerCase = WorkerCase::Piston;
        return true;
    }
    if (text == "strong-piston") {
        workerCase = WorkerCase::StrongPiston;
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
    if (text == "porous-flow") {
        workerCase = WorkerCase::PorousFlow;
        return true;
    }
    if (text == "moving-porous-flow") {
        workerCase = WorkerCase::MovingPorousFlow;
        return true;
    }
    if (text == "porous-sheet") {
        workerCase = WorkerCase::PorousSheet;
        return true;
    }
    if (text == "pressure-jump") {
        workerCase = WorkerCase::PressureJump;
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
    bool stepsRequested = false;
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
        } else if (argument == "--control-stdio") {
            options.controlStdio = true;
        } else if (argument == "--mimetic-pressure-audit") {
            options.mimeticPressureAudit = true;
        } else if (argument == "--case") {
            if (++index >= argc
                || !parseWorkerCase(argv[index], options.workerCase)) {
                error = "--case requires 'structural', 'hemisphere', 'flag', 'ram-cell', 'pressure-cell', 'piston', "
                    "'strong-piston', 'open-piston', 'periodic-flow', 'porous-flow', "
                    "'moving-porous-flow', "
                    "'porous-sheet', or "
                    "'pressure-jump'";
                return false;
            }
        } else if (argument.starts_with("--case=")) {
            if (!parseWorkerCase(argument.substr(7), options.workerCase)) {
                error = "--case requires 'structural', 'hemisphere', 'flag', 'ram-cell', 'pressure-cell', 'piston', "
                    "'strong-piston', 'open-piston', 'periodic-flow', 'porous-flow', "
                    "'moving-porous-flow', "
                    "'porous-sheet', or "
                    "'pressure-jump'";
                return false;
            }
        } else if (argument == "--steps") {
            stepsRequested = true;
            if (++index >= argc
                || !parseUnsigned(argv[index], options.steps)) {
                error = "--steps requires an unsigned integer";
                return false;
            }
        } else if (argument.starts_with("--steps=")) {
            stepsRequested = true;
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
        && options.workerCase != WorkerCase::PeriodicFlow
        && options.workerCase != WorkerCase::StrongPiston
        && options.workerCase != WorkerCase::OpenPiston
        && options.workerCase != WorkerCase::ScenePressureCell
        && options.workerCase != WorkerCase::MovingPorousFlow
        && options.workerCase != WorkerCase::PorousSheet) {
        error = "checkpoint paths require --case pressure-cell, strong-piston, "
            "moving-porous-flow, open-piston, porous-sheet, or "
            "periodic-flow";
        return false;
    }
    if (options.checkpointEvery != 0
        && options.checkpointOutputPath.empty()) {
        error = "--checkpoint-every requires --checkpoint-out";
        return false;
    }
    if (options.mimeticPressureAudit
        && options.workerCase != WorkerCase::ScenePressureCell) {
        error = "--mimetic-pressure-audit requires --case pressure-cell";
        return false;
    }
    if (options.controlStdio) {
        if (options.workerCase != WorkerCase::PeriodicFlow
            && options.workerCase != WorkerCase::StrongPiston
            && options.workerCase != WorkerCase::OpenPiston
            && options.workerCase != WorkerCase::MovingPorousFlow
            && options.workerCase != WorkerCase::PorousSheet) {
            error = "--control-stdio requires --case strong-piston, "
                "moving-porous-flow, open-piston, porous-sheet, or "
                "periodic-flow";
            return false;
        }
        if (stepsRequested) {
            error = "--control-stdio does not accept --steps";
            return false;
        }
        if (options.checkpointEvery != 0) {
            error = "--control-stdio does not accept --checkpoint-every";
            return false;
        }
        if (viewerRequested) {
            error = "--control-stdio does not launch the viewer";
            return false;
        }
        options.viewer = false;
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

bool readCheckpointBytes(const std::filesystem::path& path,
                         const std::uint64_t maximumBytes,
                         std::vector<std::uint8_t>& bytes,
                         std::string& error) {
    std::error_code fileError;
    const std::uintmax_t size = std::filesystem::file_size(path, fileError);
    if (fileError) {
        error = "cannot inspect checkpoint: " + path.string();
        return false;
    }
    if (size > maximumBytes
        || size > std::numeric_limits<std::size_t>::max()) {
        error = "checkpoint exceeds the byte limit";
        return false;
    }
    std::vector<std::uint8_t> candidate(static_cast<std::size_t>(size));
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        error = "cannot open checkpoint: " + path.string();
        return false;
    }
    if (!candidate.empty()) {
        input.read(
            reinterpret_cast<char*>(candidate.data()),
            static_cast<std::streamsize>(candidate.size()));
    }
    if (!input || input.peek() != std::char_traits<char>::eof()) {
        error = "failed to read complete checkpoint";
        return false;
    }
    bytes = std::move(candidate);
    return true;
}

bool atomicallyWriteCheckpointBytes(
    const std::filesystem::path& path,
    const std::span<const std::uint8_t> bytes,
    std::string& error) {
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
            error = "failed to write complete checkpoint";
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
        error = "cannot atomically replace checkpoint, error "
            + std::to_string(GetLastError());
        std::error_code removeError;
        std::filesystem::remove(temporary, removeError);
        return false;
    }
#else
    std::error_code renameError;
    std::filesystem::rename(temporary, path, renameError);
    if (renameError) {
        error = "cannot atomically replace checkpoint: "
            + renameError.message();
        std::error_code removeError;
        std::filesystem::remove(temporary, removeError);
        return false;
    }
#endif
    return true;
}

bool readPeriodicFlowCheckpoint(
    const std::filesystem::path& path,
    simwing::fsi::PeriodicFlowCaseCheckpoint& checkpoint,
    std::string& error) {
    const simwing::fsi::PeriodicFlowCaseCheckpointLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readCheckpointBytes(path, limits.maximumBytes, bytes, error)) {
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
    return atomicallyWriteCheckpointBytes(path, bytes, error);
}

bool readMovingPorousFlowCheckpoint(
    const std::filesystem::path& path,
    simwing::fsi::MovingPorousFlowCaseCheckpoint& checkpoint,
    std::string& error) {
    const simwing::fsi::MovingPorousFlowCaseCheckpointPersistenceLimits
        limits;
    std::vector<std::uint8_t> bytes;
    if (!readCheckpointBytes(
            path, limits.maximumEncodedBytes, bytes, error)) {
        return false;
    }
    simwing::fsi::MovingPorousFlowCaseCheckpointPersistenceError
        protocolError;
    if (!simwing::fsi::deserializeMovingPorousFlowCaseCheckpoint(
            bytes, checkpoint, &protocolError, limits)) {
        error = protocolError.message;
        return false;
    }
    return true;
}

bool writeMovingPorousFlowCheckpoint(
    const std::filesystem::path& path,
    const simwing::fsi::MovingPorousFlowCaseCheckpoint& checkpoint,
    std::string& error) {
    std::vector<std::uint8_t> bytes;
    simwing::fsi::MovingPorousFlowCaseCheckpointPersistenceError
        protocolError;
    if (!simwing::fsi::serializeMovingPorousFlowCaseCheckpoint(
            checkpoint, bytes, &protocolError)) {
        error = protocolError.message;
        return false;
    }
    return atomicallyWriteCheckpointBytes(path, bytes, error);
}

bool readOpenPistonCheckpoint(
    const std::filesystem::path& path,
    simwing::fsi::OpenPistonCaseCheckpoint& checkpoint,
    std::string& error) {
    const simwing::fsi::OpenPistonCheckpointPersistenceLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readCheckpointBytes(
            path, limits.maximumEncodedBytes, bytes, error)) {
        return false;
    }
    simwing::fsi::OpenPistonCheckpointPersistenceError protocolError;
    if (!simwing::fsi::deserializeOpenPistonCheckpoint(
            bytes, checkpoint, &protocolError, limits)) {
        error = protocolError.message;
        return false;
    }
    return true;
}

bool writeOpenPistonCheckpoint(
    const std::filesystem::path& path,
    const simwing::fsi::OpenPistonCaseCheckpoint& checkpoint,
    std::string& error) {
    std::vector<std::uint8_t> bytes;
    simwing::fsi::OpenPistonCheckpointPersistenceError protocolError;
    if (!simwing::fsi::serializeOpenPistonCheckpoint(
            checkpoint, bytes, &protocolError)) {
        error = protocolError.message;
        return false;
    }
    return atomicallyWriteCheckpointBytes(path, bytes, error);
}

bool readStrongPistonCheckpoint(
    const std::filesystem::path& path,
    simwing::fsi::StrongCoupledPistonCheckpoint& checkpoint,
    std::string& error) {
    const simwing::fsi::StrongPistonCheckpointPersistenceLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readCheckpointBytes(
            path, limits.maximumEncodedBytes, bytes, error)) {
        return false;
    }
    simwing::fsi::StrongPistonCheckpointPersistenceError protocolError;
    if (!simwing::fsi::deserializeStrongPistonCheckpoint(
            bytes, checkpoint, &protocolError, limits)) {
        error = protocolError.message;
        return false;
    }
    return true;
}

bool writeStrongPistonCheckpoint(
    const std::filesystem::path& path,
    const simwing::fsi::StrongCoupledPistonCheckpoint& checkpoint,
    std::string& error) {
    std::vector<std::uint8_t> bytes;
    simwing::fsi::StrongPistonCheckpointPersistenceError protocolError;
    if (!simwing::fsi::serializeStrongPistonCheckpoint(
            checkpoint, bytes, &protocolError)) {
        error = protocolError.message;
        return false;
    }
    return atomicallyWriteCheckpointBytes(path, bytes, error);
}

bool readScenePressureCellCheckpoint(
    const std::filesystem::path& path,
    simwing::fsi::ScenePressureCellCheckpoint& checkpoint,
    std::string& error) {
    const simwing::fsi::ScenePressureCellCheckpointPersistenceLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readCheckpointBytes(
            path, limits.maximumEncodedBytes, bytes, error)) {
        return false;
    }
    simwing::fsi::ScenePressureCellCheckpointPersistenceError protocolError;
    if (!simwing::fsi::deserializeScenePressureCellCheckpoint(
            bytes, checkpoint, &protocolError, limits)) {
        error = protocolError.message;
        return false;
    }
    return true;
}

bool writeScenePressureCellCheckpoint(
    const std::filesystem::path& path,
    const simwing::fsi::ScenePressureCellCheckpoint& checkpoint,
    std::string& error) {
    std::vector<std::uint8_t> bytes;
    simwing::fsi::ScenePressureCellCheckpointPersistenceError protocolError;
    if (!simwing::fsi::serializeScenePressureCellCheckpoint(
            checkpoint, bytes, &protocolError)) {
        error = protocolError.message;
        return false;
    }
    return atomicallyWriteCheckpointBytes(path, bytes, error);
}

bool readPorousSheetCheckpoint(
    const std::filesystem::path& path,
    simwing::fsi::CoupledPorousSheetCheckpoint& checkpoint,
    std::string& error) {
    const simwing::fsi::CoupledPorousSheetCheckpointPersistenceLimits limits;
    std::vector<std::uint8_t> bytes;
    if (!readCheckpointBytes(
            path, limits.maximumEncodedBytes, bytes, error)) {
        return false;
    }
    simwing::fsi::CoupledPorousSheetCase owner;
    simwing::fsi::CoupledPorousSheetCheckpointPersistenceError protocolError;
    if (!simwing::fsi::deserializeCoupledPorousSheetCheckpoint(
            bytes, owner, checkpoint, &protocolError, limits)) {
        error = protocolError.message;
        return false;
    }
    return true;
}

bool writePorousSheetCheckpoint(
    const std::filesystem::path& path,
    const simwing::fsi::CoupledPorousSheetCase& owner,
    const simwing::fsi::CoupledPorousSheetCheckpoint& checkpoint,
    std::string& error) {
    std::vector<std::uint8_t> bytes;
    simwing::fsi::CoupledPorousSheetCheckpointPersistenceError protocolError;
    if (!simwing::fsi::serializeCoupledPorousSheetCheckpoint(
            owner, checkpoint, bytes, &protocolError)) {
        error = protocolError.message;
        return false;
    }
    return atomicallyWriteCheckpointBytes(path, bytes, error);
}

bool configureBinaryControlStdio(std::string& error) {
#ifdef _WIN32
    if (_setmode(_fileno(stdin), _O_BINARY) == -1) {
        error = "cannot switch control stdin to binary mode";
        return false;
    }
    if (_setmode(_fileno(stdout), _O_BINARY) == -1) {
        error = "cannot switch control stdout to binary mode";
        return false;
    }
#else
    (void)error;
#endif
    return true;
}

template<typename Simulation>
int runWorkerControl(
    const Options& options,
    Simulation& simulation,
    std::ofstream& traceOutput) {
    std::string setupError;
    if (!configureBinaryControlStdio(setupError)) {
        std::fprintf(stderr, "%s\n", setupError.c_str());
        return 1;
    }

    simwing::viewer::TraceWriter traceWriter(traceOutput);
    if (!traceWriter.writeHeader(simulation.traceHeader())
        || !flushTrace(traceOutput, setupError)) {
        if (setupError.empty()) {
            setupError = traceWriter.error().message;
        }
        std::fprintf(stderr, "trace header failed: %s\n",
                     setupError.c_str());
        return 1;
    }

    bool traceFinished = false;
    const auto finishTrace = [&](std::string& finishError) {
        if (traceFinished) {
            return true;
        }
        if (!traceWriter.finish()
            || !flushTrace(traceOutput, finishError)) {
            if (finishError.empty()) {
                finishError = traceWriter.error().message;
            }
            return false;
        }
        traceFinished = true;
        return true;
    };

    using ControlHooks = std::conditional_t<
        std::is_same_v<Simulation, simwing::fsi::PeriodicFlowCase>,
        simwing::fsi::PeriodicFlowControlHooks,
        std::conditional_t<
            std::is_same_v<Simulation,
                           simwing::fsi::MovingPorousFlowCase>,
            simwing::fsi::MovingPorousFlowControlHooks,
            std::conditional_t<
                std::is_same_v<Simulation, simwing::fsi::OpenPistonCase>,
                simwing::fsi::OpenPistonControlHooks,
                std::conditional_t<
                    std::is_same_v<
                        Simulation,
                        simwing::fsi::StrongCoupledPistonWorkerCase>,
                    simwing::fsi::StrongPistonControlHooks,
                    simwing::fsi::PorousSheetControlHooks>>>>;
    ControlHooks hooks;
    hooks.publishFrame = [&](const simwing::viewer::DiagnosticFrame& frame,
                             std::string& hookError) {
        if (!traceWriter.writeFrame(frame)
            || !flushTrace(traceOutput, hookError)) {
            if (hookError.empty()) {
                hookError = traceWriter.error().message;
            }
            return false;
        }
        return true;
    };
    hooks.writeCheckpoint = [&](const auto& checkpoint,
                                std::string& hookError) {
        if (options.checkpointOutputPath.empty()) {
            hookError = "--checkpoint-out is not configured";
            return false;
        }
        if constexpr (std::is_same_v<
                          Simulation,
                          simwing::fsi::PeriodicFlowCase>) {
            return writePeriodicFlowCheckpoint(
                options.checkpointOutputPath, checkpoint, hookError);
        } else if constexpr (std::is_same_v<
                                 Simulation,
                                 simwing::fsi::MovingPorousFlowCase>) {
            return writeMovingPorousFlowCheckpoint(
                options.checkpointOutputPath, checkpoint, hookError);
        } else if constexpr (std::is_same_v<
                                 Simulation,
                                 simwing::fsi::OpenPistonCase>) {
            return writeOpenPistonCheckpoint(
                options.checkpointOutputPath, checkpoint, hookError);
        } else if constexpr (std::is_same_v<
                                 Simulation,
                                 simwing::fsi::StrongCoupledPistonWorkerCase>) {
            return writeStrongPistonCheckpoint(
                options.checkpointOutputPath, checkpoint, hookError);
        } else {
            return writePorousSheetCheckpoint(
                options.checkpointOutputPath,
                simulation,
                checkpoint,
                hookError);
        }
    };
    using ControlSession = std::conditional_t<
        std::is_same_v<Simulation, simwing::fsi::PeriodicFlowCase>,
        simwing::fsi::PeriodicFlowControlSession,
        std::conditional_t<
            std::is_same_v<Simulation,
                           simwing::fsi::MovingPorousFlowCase>,
            simwing::fsi::MovingPorousFlowControlSession,
            std::conditional_t<
                std::is_same_v<Simulation, simwing::fsi::OpenPistonCase>,
                simwing::fsi::OpenPistonControlSession,
                std::conditional_t<
                    std::is_same_v<
                        Simulation,
                        simwing::fsi::StrongCoupledPistonWorkerCase>,
                    simwing::fsi::StrongPistonControlSession,
                    simwing::fsi::PorousSheetControlSession>>>>;
    ControlSession session(
        simulation, std::move(hooks));

    simwing::fsi::WorkerControlStreamError streamError;
    if (!simwing::fsi::writeWorkerControlResponse(
            std::cout, session.readyResponse(), &streamError)) {
        std::string finishError;
        (void)finishTrace(finishError);
        std::fprintf(stderr, "control response failed: %s\n",
                     streamError.message.c_str());
        return 1;
    }

    for (;;) {
        simwing::fsi::WorkerControlCommand command;
        const simwing::fsi::WorkerControlStreamResult readResult =
            simwing::fsi::readWorkerControlCommand(
                std::cin, command, &streamError);
        if (readResult
            == simwing::fsi::WorkerControlStreamResult::EndOfStream) {
            std::string finishError;
            if (!finishTrace(finishError)) {
                std::fprintf(stderr, "trace completion failed: %s\n",
                             finishError.c_str());
            }
            std::fprintf(
                stderr,
                "control stream ended before an explicit stop command\n");
            return 1;
        }
        if (readResult == simwing::fsi::WorkerControlStreamResult::Error) {
            std::string finishError;
            if (!finishTrace(finishError)) {
                std::fprintf(stderr, "trace completion failed: %s\n",
                             finishError.c_str());
            }
            std::fprintf(stderr, "control command failed: %s\n",
                         streamError.message.c_str());
            return 2;
        }

        simwing::fsi::WorkerControlResponse response;
        simwing::fsi::WorkerControlProtocolError protocolError;
        if (!session.execute(command, response, &protocolError)) {
            std::string finishError;
            if (!finishTrace(finishError)) {
                std::fprintf(stderr, "trace completion failed: %s\n",
                             finishError.c_str());
            }
            std::fprintf(stderr, "control command is invalid: %s\n",
                         protocolError.message.c_str());
            return 2;
        }

        int terminalResult = 0;
        if (response.kind
            == simwing::fsi::WorkerControlResponseKind::Stopped) {
            std::string finishError;
            if (!finishTrace(finishError)) {
                response.kind =
                    simwing::fsi::WorkerControlResponseKind::Error;
                response.failureCode =
                    simwing::fsi::WorkerControlFailureCode::InternalFailure;
                response.errorMessage = "trace completion failed";
                if (!finishError.empty()) {
                    response.errorMessage += ": " + finishError;
                }
                const std::size_t limit =
                    simwing::fsi::WorkerControlProtocolLimits{}
                        .maximumErrorMessageBytes;
                if (response.errorMessage.size() > limit) {
                    response.errorMessage.resize(limit);
                }
                terminalResult = 1;
            }
        }

        if (!simwing::fsi::writeWorkerControlResponse(
                std::cout, response, &streamError)) {
            std::string finishError;
            (void)finishTrace(finishError);
            std::fprintf(stderr, "control response failed: %s\n",
                         streamError.message.c_str());
            return 1;
        }
        if (session.stopped()) {
            return terminalResult;
        }
    }
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
        std::optional<simwing::fsi::OpenPistonCaseCheckpoint>
            restoredOpenPistonCheckpoint;
        std::optional<simwing::fsi::StrongCoupledPistonCheckpoint>
            restoredStrongPistonCheckpoint;
        std::optional<simwing::fsi::ScenePressureCellCheckpoint>
            restoredScenePressureCellCheckpoint;
        std::optional<simwing::fsi::MovingPorousFlowCaseCheckpoint>
            restoredMovingPorousFlowCheckpoint;
        std::optional<simwing::fsi::CoupledPorousSheetCheckpoint>
            restoredPorousSheetCheckpoint;
        if (!options.checkpointInputPath.empty()) {
            bool restored = false;
            if (options.workerCase == WorkerCase::PeriodicFlow) {
                restoredPeriodicFlowCheckpoint.emplace();
                restored = readPeriodicFlowCheckpoint(
                    options.checkpointInputPath,
                    *restoredPeriodicFlowCheckpoint, error);
            } else if (options.workerCase == WorkerCase::OpenPiston) {
                restoredOpenPistonCheckpoint.emplace();
                restored = readOpenPistonCheckpoint(
                    options.checkpointInputPath,
                    *restoredOpenPistonCheckpoint, error);
            } else if (options.workerCase == WorkerCase::StrongPiston) {
                restoredStrongPistonCheckpoint.emplace();
                restored = readStrongPistonCheckpoint(
                    options.checkpointInputPath,
                    *restoredStrongPistonCheckpoint, error);
            } else if (options.workerCase
                       == WorkerCase::ScenePressureCell) {
                restoredScenePressureCellCheckpoint.emplace();
                restored = readScenePressureCellCheckpoint(
                    options.checkpointInputPath,
                    *restoredScenePressureCellCheckpoint, error);
            } else if (options.workerCase
                       == WorkerCase::MovingPorousFlow) {
                restoredMovingPorousFlowCheckpoint.emplace();
                restored = readMovingPorousFlowCheckpoint(
                    options.checkpointInputPath,
                    *restoredMovingPorousFlowCheckpoint, error);
            } else {
                restoredPorousSheetCheckpoint.emplace();
                restored = readPorousSheetCheckpoint(
                    options.checkpointInputPath,
                    *restoredPorousSheetCheckpoint, error);
            }
            if (!restored) {
                std::fprintf(stderr, "checkpoint restore failed: %s\n",
                             error.c_str());
                return 1;
            }
        }

        if (!options.controlStdio
            && options.workerCase == WorkerCase::MovingPorousFlow
            && !options.checkpointOutputPath.empty()) {
            const simwing::fsi::MovingPorousFlowCaseCheckpointPersistenceLimits
                limits;
            const std::uint64_t restoredSteps =
                restoredMovingPorousFlowCheckpoint
                ? restoredMovingPorousFlowCheckpoint->acceptedStepCount
                : 0;
            if (restoredSteps > limits.maximumReplaySteps
                || options.steps
                    > limits.maximumReplaySteps - restoredSteps) {
                std::fprintf(
                    stderr,
                    "moving-porous-flow checkpoint output would exceed "
                    "the %llu-step deterministic replay limit\n",
                    static_cast<unsigned long long>(
                        limits.maximumReplaySteps));
                return 2;
            }
        }

        std::ofstream output(options.tracePath,
                             std::ios::binary | std::ios::trunc);
        if (!output) {
            std::fprintf(stderr, "cannot open trace for writing: %s\n",
                         options.tracePath.string().c_str());
            return 1;
        }

        if (options.controlStdio) {
            if (options.workerCase == WorkerCase::PeriodicFlow) {
                simwing::fsi::PeriodicFlowCase simulation;
                if (restoredPeriodicFlowCheckpoint) {
                    simulation.restore(*restoredPeriodicFlowCheckpoint);
                }
                return runWorkerControl(options, simulation, output);
            }
            if (options.workerCase == WorkerCase::OpenPiston) {
                simwing::fsi::OpenPistonCase simulation;
                if (restoredOpenPistonCheckpoint) {
                    simulation.restore(*restoredOpenPistonCheckpoint);
                }
                return runWorkerControl(options, simulation, output);
            }
            if (options.workerCase == WorkerCase::StrongPiston) {
                simwing::fsi::StrongCoupledPistonWorkerCase simulation;
                if (restoredStrongPistonCheckpoint) {
                    simulation.restore(*restoredStrongPistonCheckpoint);
                }
                return runWorkerControl(options, simulation, output);
            }
            if (options.workerCase == WorkerCase::MovingPorousFlow) {
                simwing::fsi::MovingPorousFlowCase simulation;
                if (restoredMovingPorousFlowCheckpoint) {
                    simulation.restore(
                        *restoredMovingPorousFlowCheckpoint);
                }
                return runWorkerControl(options, simulation, output);
            }
            simwing::fsi::CoupledPorousSheetCase simulation;
            if (restoredPorousSheetCheckpoint) {
                simulation.restore(*restoredPorousSheetCheckpoint);
            }
            return runWorkerControl(options, simulation, output);
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
                using Simulation =
                    std::remove_cvref_t<decltype(simulation)>;
                if constexpr (
                    std::is_same_v<Simulation,
                                   simwing::fsi::PeriodicFlowCase>
                    || std::is_same_v<Simulation,
                                       simwing::fsi::OpenPistonCase>
                    || std::is_same_v<Simulation,
                                      simwing::fsi::StrongCoupledPistonWorkerCase>
                    || std::is_same_v<Simulation,
                                      simwing::fsi::ScenePressureCellCase>
                    || std::is_same_v<Simulation,
                                      simwing::fsi::MovingPorousFlowCase>
                    || std::is_same_v<Simulation,
                                      simwing::fsi::CoupledPorousSheetCase>) {
                    if (options.checkpointEvery != 0
                        && simulation.acceptedStepCount()
                               % options.checkpointEvery == 0) {
                        const bool saved = [&] {
                            if constexpr (std::is_same_v<
                                              Simulation,
                                              simwing::fsi::PeriodicFlowCase>) {
                                return writePeriodicFlowCheckpoint(
                                    options.checkpointOutputPath,
                                    simulation.checkpoint(), error);
                            } else if constexpr (std::is_same_v<
                                                     Simulation,
                                                     simwing::fsi::OpenPistonCase>) {
                                return writeOpenPistonCheckpoint(
                                    options.checkpointOutputPath,
                                    simulation.checkpoint(), error);
                            } else if constexpr (std::is_same_v<
                                                     Simulation,
                                                     simwing::fsi::StrongCoupledPistonWorkerCase>) {
                                return writeStrongPistonCheckpoint(
                                    options.checkpointOutputPath,
                                    simulation.checkpoint(), error);
                            } else if constexpr (std::is_same_v<
                                                     Simulation,
                                                     simwing::fsi::ScenePressureCellCase>) {
                                return writeScenePressureCellCheckpoint(
                                    options.checkpointOutputPath,
                                    simulation.checkpoint(), error);
                            } else if constexpr (std::is_same_v<
                                                     Simulation,
                                                     simwing::fsi::MovingPorousFlowCase>) {
                                return writeMovingPorousFlowCheckpoint(
                                    options.checkpointOutputPath,
                                    simulation.checkpoint(), error);
                            } else {
                                return writePorousSheetCheckpoint(
                                    options.checkpointOutputPath,
                                    simulation,
                                    simulation.checkpoint(), error);
                            }
                        }();
                        if (!saved) {
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

            using Simulation =
                std::remove_cvref_t<decltype(simulation)>;
            if constexpr (
                std::is_same_v<Simulation,
                               simwing::fsi::PeriodicFlowCase>
                || std::is_same_v<Simulation,
                                   simwing::fsi::OpenPistonCase>
                || std::is_same_v<Simulation,
                                  simwing::fsi::StrongCoupledPistonWorkerCase>
                || std::is_same_v<Simulation,
                                  simwing::fsi::ScenePressureCellCase>
                || std::is_same_v<Simulation,
                                  simwing::fsi::MovingPorousFlowCase>
                || std::is_same_v<Simulation,
                                  simwing::fsi::CoupledPorousSheetCase>) {
                if (!options.checkpointOutputPath.empty()
                    && lastCheckpointStep
                        != simulation.acceptedStepCount()) {
                    const bool saved = [&] {
                        if constexpr (std::is_same_v<
                                          Simulation,
                                          simwing::fsi::PeriodicFlowCase>) {
                            return writePeriodicFlowCheckpoint(
                                options.checkpointOutputPath,
                                simulation.checkpoint(), error);
                        } else if constexpr (std::is_same_v<
                                                 Simulation,
                                                 simwing::fsi::OpenPistonCase>) {
                            return writeOpenPistonCheckpoint(
                                options.checkpointOutputPath,
                                simulation.checkpoint(), error);
                        } else if constexpr (std::is_same_v<
                                                 Simulation,
                                                 simwing::fsi::StrongCoupledPistonWorkerCase>) {
                            return writeStrongPistonCheckpoint(
                                options.checkpointOutputPath,
                                simulation.checkpoint(), error);
                        } else if constexpr (std::is_same_v<
                                                 Simulation,
                                                 simwing::fsi::ScenePressureCellCase>) {
                            return writeScenePressureCellCheckpoint(
                                options.checkpointOutputPath,
                                simulation.checkpoint(), error);
                        } else if constexpr (std::is_same_v<
                                                 Simulation,
                                                 simwing::fsi::MovingPorousFlowCase>) {
                            return writeMovingPorousFlowCheckpoint(
                                options.checkpointOutputPath,
                                simulation.checkpoint(), error);
                        } else {
                            return writePorousSheetCheckpoint(
                                options.checkpointOutputPath,
                                simulation,
                                simulation.checkpoint(), error);
                        }
                    }();
                    if (!saved) {
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
                if constexpr (std::is_same_v<
                                  Simulation,
                                  simwing::fsi::StrongCoupledPistonWorkerCase>) {
                    const auto& coupled = simulation.diagnostics();
                    std::printf(
                        "simwing-fsi completed %llu strong-piston step(s), "
                        "t=%.9g s, speed=%.9g m/s, coupling-iterations=%llu, "
                        "solver-runs=%llu, attempts=%zu, retries=%u, "
                        "velocity-closure=%.3g m/s, "
                        "added-mass=%.9g kg, analytic-residual=%.3g m/s, "
                        "checkpoint-writes=%llu, "
                        "trace=%s\n",
                        static_cast<unsigned long long>(
                            checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        coupled.acceptedSpeedMetersPerSecond,
                        static_cast<unsigned long long>(
                            coupled.coupling.lastIteration
                                .convergence.iteration),
                        static_cast<unsigned long long>(
                            coupled.coupling.solverRunCount),
                        coupled.coupling.attempts.size(),
                        coupled.coupling.decision.retryCount,
                        coupled.velocityClosureMetersPerSecond,
                        coupled.measuredDiscreteAddedMassKilograms,
                        coupled.analyticSpeedResidualMetersPerSecond,
                        static_cast<unsigned long long>(
                            checkpointWriteCount),
                        options.tracePath.string().c_str());
                } else if constexpr (std::is_same_v<
                                  Simulation,
                                  simwing::fsi::OpenPistonCase>) {
                    std::printf(
                        "simwing-fsi completed %llu open-piston step(s), "
                        "t=%.9g s, "
                        "bounds=[%.6g %.6g %.6g]-[%.6g %.6g %.6g] m, "
                        "max-strain=%.6g, checkpoint-writes=%llu, "
                        "trace=%s\n",
                        static_cast<unsigned long long>(
                            checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        minimum.x, minimum.y, minimum.z,
                        maximum.x, maximum.y, maximum.z,
                        diagnostics.maximumAbsoluteMembraneStrain,
                        static_cast<unsigned long long>(
                            checkpointWriteCount),
                        options.tracePath.string().c_str());
                } else if constexpr (std::is_same_v<
                                         Simulation,
                                         simwing::fsi::CoupledPorousSheetCase>) {
                    const auto& coupled = simulation.diagnostics();
                    const double sheetPosition =
                        checkpoint.nodes.front().positionMeters.x;
                    const double sheetSpeed =
                        checkpoint.nodes.front().velocityMetersPerSecond.x;
                    const double fluidSpeed =
                        simulation.velocity().xFaces().front();
                    std::printf(
                        "simwing-fsi completed %llu porous-sheet step(s), "
                        "t=%.9g s, sheet-x=%.9g m, sheet-speed=%.9g m/s, "
                        "fluid-speed=%.9g m/s, energy-residual=%.3g J, "
                        "topology-rebases=%llu, face=%zu, image=%lld, "
                        "checkpoint-writes=%llu, trace=%s\n",
                        static_cast<unsigned long long>(
                        checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        sheetPosition,
                        sheetSpeed,
                        fluidSpeed,
                        coupled.energyResidualJoules,
                        static_cast<unsigned long long>(
                            simulation.topologyRebaseCount()),
                        simulation.porousFaceCoordinate(),
                        static_cast<long long>(
                            simulation.porousTopology().periodicImage),
                        static_cast<unsigned long long>(
                            checkpointWriteCount),
                        options.tracePath.string().c_str());
                } else if constexpr (std::is_same_v<
                                         Simulation,
                                         simwing::fsi::ProjectedGustFlagCase>) {
                    const auto& transfer =
                        simulation.transferDiagnostics();
                    std::printf(
                        "simwing-fsi completed %llu flag step(s), "
                        "t=%.9g s, gust=%.6g m/s, normal-motion=%.6g m, "
                        "free-edge-motion=%.6g m, pressure-force-x=%.6g N, "
                        "direct-force-x=%.6g N, reaction-force-x=%.6g N, "
                        "transfer-residual=%.3g N, trace=%s\n",
                        static_cast<unsigned long long>(
                            checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        simulation.gustSpeedMetersPerSecond(),
                        simulation.maximumNormalDisplacementMeters(),
                        simulation.maximumFreeEdgeDisplacementMeters(),
                        transfer.fluidPressureForceNewtons.x,
                        transfer.fluidLoadForceNewtons.x
                            - transfer.fluidPressureForceNewtons.x,
                        transfer.fluidLoadForceNewtons.x,
                        transfer.forceResidualNormNewtons,
                        options.tracePath.string().c_str());
                } else if constexpr (std::is_same_v<
                                         Simulation,
                                         simwing::fsi::RamAirCellCase>) {
                    const auto& coupled = simulation.diagnostics();
                    std::printf(
                        "simwing-fsi completed %llu ram-cell step(s), "
                        "t=%.9g s, gust=%.6g m/s, mouth-mean=%.3g m/s, "
                        "mouth-rms=%.6g m/s, max-motion=%.6g m, "
                        "outward-inflation=%.6g m, "
                        "reaction=[%.6g %.6g %.6g] N, "
                        "transfer-residual=%.3g N, divergence=%.3g 1/s, "
                        "trace=%s\n",
                        static_cast<unsigned long long>(
                            checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        simulation.gustSpeedMetersPerSecond(),
                        simulation.openingMeanVelocityMetersPerSecond(),
                        simulation.openingRmsVelocityMetersPerSecond(),
                        simulation.maximumDisplacementMeters(),
                        simulation.maximumOutwardInflationMeters(),
                        coupled.fluidReactionForceNewtons.x,
                        coupled.fluidReactionForceNewtons.y,
                        coupled.fluidReactionForceNewtons.z,
                        std::sqrt(
                            coupled.forceResidualNewtons.x
                                * coupled.forceResidualNewtons.x
                            + coupled.forceResidualNewtons.y
                                * coupled.forceResidualNewtons.y
                            + coupled.forceResidualNewtons.z
                                * coupled.forceResidualNewtons.z),
                        coupled.fluidDivergenceL2PerSecond,
                        options.tracePath.string().c_str());
                } else if constexpr (std::is_same_v<
                                         Simulation,
                                         simwing::fsi::ScenePressureCellCase>) {
                    const auto& coupled = simulation.diagnostics();
                    std::printf(
                        "simwing-fsi completed %llu pressure-cell step(s), "
                        "t=%.9g s, wind=%.6g m/s, pump=%.6g N, pressure-force="
                        "[%.6g %.6g %.6g] N, wall-force="
                        "[%.6g %.6g %.6g] N, max-pressure=%.6g Pa, "
                        "max-motion=%.6g m, coupling-iterations=%llu, "
                        "mac-speed=%.6g m/s, subface-spread=%.3g m/s, "
                        "bulk-viscous-loss=%.3g J, region-loss=%.3g J, "
                        "region-gcl=%.3g m^3, region-momentum-residual="
                        "%.3g N s, wall-loss=%.3g J, "
                        "wall-momentum-residual=%.3g N s, "
                        "load-closure=%.3g N, "
                        "checkpoint-writes=%llu, "
                        "trace=%s\n",
                        static_cast<unsigned long long>(
                            checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        coupled.targetMeanWindMetersPerSecond,
                        coupled.flowPumpForceNewtons,
                        coupled.pressureForceNewtons.x,
                        coupled.pressureForceNewtons.y,
                        coupled.pressureForceNewtons.z,
                        coupled.wallForceNewtons.x,
                        coupled.wallForceNewtons.y,
                        coupled.wallForceNewtons.z,
                        coupled.maximumAbsolutePressurePascals,
                        coupled.maximumDisplacementMeters,
                        static_cast<unsigned long long>(
                            coupled.coupling.solverRunCount),
                        coupled.macVelocity
                            .maximumAbsoluteVelocityMetersPerSecond,
                        coupled.macVelocity
                            .maximumSubfaceVelocityDeviationMetersPerSecond,
                        coupled.bulkFlow.firstHalfViscousEnergyLossJoules
                            + coupled.bulkFlow
                                .secondHalfViscousEnergyLossJoules,
                        coupled.regionTransport.advectiveEnergyLossJoules
                            + coupled.regionTransport.viscousEnergyLossJoules,
                        coupled.regionTransport
                            .maximumAbsoluteGeometryVolumeChangeCubicMeters,
                        coupled.regionTransport
                            .momentumResidualNormKilogramMetersPerSecond,
                        coupled.coupling.regionWall
                            .viscousDissipationJoules,
                        coupled.coupling.regionWall
                            .momentumResidualNormKilogramMetersPerSecond,
                        coupled.coupling.interfaceForceClosureNewtons,
                        static_cast<unsigned long long>(
                            checkpointWriteCount),
                        options.tracePath.string().c_str());
                    if (options.mimeticPressureAudit) {
                        const auto* audit =
                            simulation.acceptedMimeticPressureAudit();
                        const auto* comparison =
                            simulation.acceptedMimeticPressureComparison();
                        if (audit == nullptr || comparison == nullptr) {
                            std::printf(
                                "simwing-fsi mimetic-pressure-audit not run\n");
                        } else {
                            std::printf(
                                "simwing-fsi mimetic-pressure-audit accepted, "
                                "controls=%zu, shared-traces=%zu, "
                                "wall-traces=%zu, iterations=%zu, "
                                "pressure-rms-delta=%.6g Pa, "
                                "pressure-relative-delta=%.6g, "
                                "force-delta=%.6g N, "
                                "force-relative-delta=%.6g, "
                                "consecutive=%u, wall-predictor=%u\n",
                                audit->controlCells.controlCells.size(),
                                audit->condensedTraceSystem.traces.size(),
                                audit->condensedTraceSystem
                                    .eliminatedMaterialWallTraceCount,
                                audit->pressureEpoch.diagnostics.pressureSolve
                                    .reducedTraceSolve.iterationCount,
                                comparison->diagnostics
                                    .pressureDifferenceDeltaRmsPascals,
                                comparison->diagnostics
                                    .relativePressureDifferenceDeltaL2,
                                comparison->diagnostics
                                    .forceDeltaNormNewtons,
                                comparison->diagnostics.relativeForceDelta,
                                audit->usesConsecutiveWarmStart ? 1U : 0U,
                                audit->usesRegionWallPrediction ? 1U : 0U);
                        }
                    }
                } else if constexpr (std::is_same_v<
                                         Simulation,
                                         simwing::fsi::AnchoredHemisphereCase>) {
                    std::printf(
                        "simwing-fsi completed %llu hemisphere step(s), "
                        "t=%.9g s, pressure=%.6g Pa, apex-motion=%.6g m, "
                        "free-rim-motion=%.6g m, "
                        "bounds=[%.6g %.6g %.6g]-[%.6g %.6g %.6g] m, "
                        "max-strain=%.6g, trace=%s\n",
                        static_cast<unsigned long long>(
                            checkpoint.acceptedStepCount),
                        checkpoint.simulationTimeSeconds,
                        simulation.pressurePascals(),
                        simulation.apexRadialDisplacementMeters(),
                        simulation.maximumFreeRimDisplacementMeters(),
                        minimum.x, minimum.y, minimum.z,
                        maximum.x, maximum.y, maximum.z,
                        diagnostics.maximumAbsoluteMembraneStrain,
                        options.tracePath.string().c_str());
                } else {
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
                }
            } else if constexpr (std::is_same_v<
                                     Simulation,
                                     simwing::fsi::PressureJumpCase>) {
                const auto& diagnostics = simulation.diagnostics();
                std::printf(
                    "simwing-fsi completed %llu pressure-jump step(s), "
                    "t=%.9g s, crossings=%zu, divergence-L2=%.6g 1/s, "
                    "energy=%.9g J, trace=%s\n",
                    static_cast<unsigned long long>(
                        simulation.acceptedStepCount()),
                    simulation.simulationTimeSeconds(),
                    simulation.pressureJumps().faceCount(),
                    diagnostics.divergenceL2AfterPerSecond,
                    diagnostics.kineticEnergyAfterJoules,
                    options.tracePath.string().c_str());
            } else if constexpr (std::is_same_v<
                                     Simulation,
                                     simwing::fsi::MovingPorousFlowCase>) {
                const auto& diagnostics = simulation.diagnostics();
                std::printf(
                    "simwing-fsi completed %llu moving-porous-flow step(s), "
                    "t=%.9g s, topology-rebases=%llu, face=%zu, image=%lld, "
                    "sheet-position=%.9g m, kinematic-residual=%.3g m, "
                    "momentum-residual=%.3g N*s, checkpoint-writes=%llu, "
                    "trace=%s\n",
                    static_cast<unsigned long long>(
                        simulation.acceptedStepCount()),
                    simulation.simulationTimeSeconds(),
                    static_cast<unsigned long long>(
                        simulation.topologyRebaseCount()),
                    simulation.porousTopology().faceCoordinate,
                    static_cast<long long>(
                        simulation.porousTopology().periodicImage),
                    simulation.sheetPositionMeters(),
                    diagnostics.kinematicResidualMeters,
                    diagnostics.flow.momentumResidualNormNewtonSeconds,
                    static_cast<unsigned long long>(checkpointWriteCount),
                    options.tracePath.string().c_str());
            } else if constexpr (std::is_same_v<
                                     Simulation,
                                     simwing::fsi::PorousFlowCase>) {
                const auto& plug = simulation.plugDiagnostics();
                std::printf(
                    "simwing-fsi completed %llu porous-flow step(s), "
                    "t=%.9g s, speed=%.9g m/s, pressure-drop=%.9g Pa, "
                    "dissipation=%.9g J, energy-residual=%.3g J, "
                    "trace=%s\n",
                    static_cast<unsigned long long>(
                        simulation.acceptedStepCount()),
                    simulation.simulationTimeSeconds(),
                    plug.velocityAfterMetersPerSecond,
                    plug.endpointPressureDropPascals,
                    plug.porousDissipationJoules,
                    plug.energyResidualJoules,
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
        if (options.workerCase == WorkerCase::StrongPiston) {
            simwing::fsi::StrongCoupledPistonWorkerCase simulation;
            if (restoredStrongPistonCheckpoint) {
                simulation.restore(*restoredStrongPistonCheckpoint);
            }
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::OpenPiston) {
            simwing::fsi::OpenPistonCase simulation;
            if (restoredOpenPistonCheckpoint) {
                simulation.restore(*restoredOpenPistonCheckpoint);
            }
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::PeriodicFlow) {
            simwing::fsi::PeriodicFlowCase simulation;
            if (restoredPeriodicFlowCheckpoint) {
                simulation.restore(*restoredPeriodicFlowCheckpoint);
            }
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::PressureJump) {
            simwing::fsi::PressureJumpCase simulation;
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::PorousFlow) {
            simwing::fsi::PorousFlowCase simulation;
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::MovingPorousFlow) {
            simwing::fsi::MovingPorousFlowCase simulation;
            if (restoredMovingPorousFlowCheckpoint) {
                simulation.restore(*restoredMovingPorousFlowCheckpoint);
            }
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::PorousSheet) {
            simwing::fsi::CoupledPorousSheetCase simulation;
            if (restoredPorousSheetCheckpoint) {
                simulation.restore(*restoredPorousSheetCheckpoint);
            }
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::Hemisphere) {
            simwing::fsi::AnchoredHemisphereCase simulation;
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::ProjectedFlag) {
            simwing::fsi::ProjectedGustFlagCase simulation;
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::RamAirCell) {
            simwing::fsi::RamAirCellCase simulation;
            return run(simulation);
        }
        if (options.workerCase == WorkerCase::ScenePressureCell) {
            simwing::fsi::ScenePressureCellCase simulation(
                options.mimeticPressureAudit);
            if (restoredScenePressureCellCheckpoint) {
                simulation.restore(*restoredScenePressureCellCheckpoint);
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
