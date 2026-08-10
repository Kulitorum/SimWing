#include "scene_pressure_cell_mimetic_conductance_phase_refinement_audit.h"

#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

constexpr std::size_t maximumResolution = 64;

struct Options {
    std::size_t resolution = 0;
    std::optional<std::size_t> phaseIndex;
    bool allPhases = false;
    bool help = false;
};

void printUsage(FILE* stream) {
    std::fprintf(
        stream,
        "Usage: simwing-mimetic-conductance-audit --resolution N "
        "(--phase N | --all-phases)\n"
        "\n"
        "Runs the immutable graph-free pressure-cell terminal-conductance "
        "audit with\n"
        "the fingerprinted offline numeric profile. Resolution is bounded "
        "to 2..64;\n"
        "phase is one canonical half-cell translation index in 0..7. "
        "--all-phases\n"
        "is deliberately explicit because fine-grid samples can take "
        "minutes. The tool\n"
        "does not construct the graph pressure operator or alter worker "
        "arithmetic.\n");
}

bool parseUnsigned(const std::string_view text, std::size_t& value) {
    if (text.empty()) {
        return false;
    }
    std::uint64_t parsed = 0;
    const char* begin = text.data();
    const char* end = begin + text.size();
    const auto [position, error] = std::from_chars(begin, end, parsed);
    if (error != std::errc{} || position != end
        || parsed > std::numeric_limits<std::size_t>::max()) {
        return false;
    }
    value = static_cast<std::size_t>(parsed);
    return true;
}

bool setResolution(const std::string_view text,
                   Options& options,
                   std::string& error) {
    std::size_t resolution = 0;
    if (options.resolution != 0 || !parseUnsigned(text, resolution)
        || resolution < 2 || resolution > maximumResolution) {
        error = "--resolution requires one integer in 2..64";
        return false;
    }
    options.resolution = resolution;
    return true;
}

bool setPhase(const std::string_view text,
              Options& options,
              std::string& error) {
    std::size_t phaseIndex = 0;
    if (options.phaseIndex.has_value() || options.allPhases
        || !parseUnsigned(text, phaseIndex)
        || phaseIndex
            >= simwing::fsi::
                scenePressureCellMimeticConductanceCanonicalGridPhases.size()) {
        error = "--phase requires one canonical phase index in 0..7 and conflicts with --all-phases";
        return false;
    }
    options.phaseIndex = phaseIndex;
    return true;
}

bool parseOptions(const int argc,
                  char* argv[],
                  Options& options,
                  std::string& error) {
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            options.help = true;
        } else if (argument == "--resolution") {
            if (++index >= argc
                || !setResolution(argv[index], options, error)) {
                if (error.empty()) {
                    error = "--resolution requires one integer in 2..64";
                }
                return false;
            }
        } else if (argument.starts_with("--resolution=")) {
            if (!setResolution(argument.substr(13), options, error)) {
                return false;
            }
        } else if (argument == "--phase") {
            if (++index >= argc || !setPhase(argv[index], options, error)) {
                if (error.empty()) {
                    error = "--phase requires one canonical phase index in 0..7";
                }
                return false;
            }
        } else if (argument.starts_with("--phase=")) {
            if (!setPhase(argument.substr(8), options, error)) {
                return false;
            }
        } else if (argument == "--all-phases") {
            if (options.allPhases || options.phaseIndex.has_value()) {
                error = "--all-phases may appear once and conflicts with --phase";
                return false;
            }
            options.allPhases = true;
        } else {
            error = "unknown option: " + std::string(argument);
            return false;
        }
    }
    if (options.help) {
        return true;
    }
    if (options.resolution == 0) {
        error = "--resolution is required";
        return false;
    }
    if (!options.phaseIndex.has_value() && !options.allPhases) {
        error = "exactly one of --phase or --all-phases is required";
        return false;
    }
    return true;
}

const char* statusName(
    const simwing::fsi::
        ScenePressureCellMimeticConductancePhaseSampleStatus status) {
    using Status = simwing::fsi::
        ScenePressureCellMimeticConductancePhaseSampleStatus;
    switch (status) {
    case Status::Accepted:
        return "accepted";
    case Status::RejectedLocalCellLinearConsistency:
        return "rejected-local-cell-linear-consistency";
    }
    return "unknown";
}

void printAudit(
    const simwing::fsi::
        ScenePressureCellMimeticConductancePhaseRefinementAudit& audit,
    const std::span<const std::size_t> canonicalPhaseIndices) {
    const auto& level = audit.levels.front();
    std::cout.imbue(std::locale::classic());
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "audit version=" << audit.version
              << " resolution=" << level.cellCounts.x
              << " phase-count=" << level.samples.size()
              << " fingerprint=0x" << std::hex << audit.fingerprint
              << " structure-fingerprint=0x"
              << audit.structureDefinitionFingerprint << std::dec << '\n';
    for (std::size_t index = 0; index < level.samples.size(); ++index) {
        const auto& sample = level.samples[index];
        const auto phase = sample.gridPhaseFraction;
        std::cout << "sample phase=" << canonicalPhaseIndices[index]
                  << " translation=(" << phase.x << ',' << phase.y << ','
                  << phase.z << ')'
                  << " status=" << statusName(sample.status)
                  << " grid-cells=" << sample.gridCellCount
                  << " controls=" << sample.controlVolumeCount
                  << " full-traces=" << sample.fullTraceCount
                  << " reduced-traces=" << sample.reducedTraceCount
                  << " opening-traces=" << sample.openingTraceCount
                  << " intake-area-m2=" << sample.intakeAreaSquareMeters
                  << " conductance-m=" << sample.conductanceMeters
                  << " normalized=" << sample.normalizedConductance;
        if (sample.conductanceAudit.has_value()) {
            const auto& solve = sample.conductanceAudit->solveDiagnostics
                                    .reducedTraceSolve;
            std::cout << " iterations=" << solve.iterationCount
                      << " residual-l2-pa-m="
                      << solve.finalResidualL2PascalsMeters;
        } else if (sample.localCellLinearConsistencyRejection.has_value()) {
            const auto& rejection =
                *sample.localCellLinearConsistencyRejection;
            std::cout << " rejected-control=" << rejection.controlCellIndex
                      << " algebraic-error="
                      << rejection.localCell
                             .maximumAlgebraicConsistencyError
                      << " algebraic-tolerance="
                      << rejection.localCell.algebraicConsistencyTolerance;
        }
        std::cout << '\n';
    }
    std::cout << "summary accepted=" << level.acceptedSampleCount
              << " rejected="
              << level.rejectedLocalCellLinearConsistencySampleCount
              << " minimum=" << level.minimumNormalizedConductance
              << " maximum=" << level.maximumNormalizedConductance
              << " mean=" << level.meanNormalizedConductance
              << " cv="
              << level.normalizedConductanceCoefficientOfVariation << '\n';
}

} // namespace

int main(const int argc, char* argv[]) {
    Options options;
    std::string error;
    if (!parseOptions(argc, argv, options, error)) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        printUsage(stderr);
        return 2;
    }
    if (options.help) {
        printUsage(stdout);
        return 0;
    }

    try {
        std::vector<std::size_t> canonicalPhaseIndices;
        std::vector<simwing::fsi::fluid::Vector3> phases;
        if (options.allPhases) {
            canonicalPhaseIndices.reserve(
                simwing::fsi::
                    scenePressureCellMimeticConductanceCanonicalGridPhases
                        .size());
            phases.reserve(canonicalPhaseIndices.capacity());
            for (std::size_t index = 0;
                 index < simwing::fsi::
                     scenePressureCellMimeticConductanceCanonicalGridPhases
                         .size();
                 ++index) {
                canonicalPhaseIndices.push_back(index);
                phases.push_back(simwing::fsi::
                    scenePressureCellMimeticConductanceCanonicalGridPhases[
                        index]);
            }
        } else {
            canonicalPhaseIndices.push_back(*options.phaseIndex);
            phases.push_back(simwing::fsi::
                scenePressureCellMimeticConductanceCanonicalGridPhases[
                    *options.phaseIndex]);
        }
        const std::array<simwing::fsi::fluid::GridCellCounts, 1>
            resolutions{{{
                options.resolution, options.resolution,
                options.resolution,
            }}};
        const auto audit = simwing::fsi::
            auditScenePressureCellMimeticConductancePhaseRefinement(
                resolutions, phases,
                simwing::fsi::
                    makeScenePressureCellMimeticConductanceOfflineAuditSettings());
        printAudit(audit, canonicalPhaseIndices);
        return 0;
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "error: audit failed: %s\n", exception.what());
        return 1;
    }
}
