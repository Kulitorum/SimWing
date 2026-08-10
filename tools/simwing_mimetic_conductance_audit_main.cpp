#include "scene_pressure_cell_mimetic_conductance_phase_refinement_audit.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <optional>
#include <stdexcept>
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

class CompensatedSum final {
public:
    void add(const double value) noexcept {
        const double next = sum_ + value;
        if (std::abs(sum_) >= std::abs(value)) {
            correction_ += (sum_ - next) + value;
        } else {
            correction_ += (value - next) + sum_;
        }
        sum_ = next;
    }

    [[nodiscard]] double value() const noexcept {
        return sum_ + correction_;
    }

private:
    double sum_ = 0.0;
    double correction_ = 0.0;
};

class CompactAggregate final {
public:
    void add(const simwing::fsi::
                 ScenePressureCellMimeticConductancePhaseSample& sample) {
        using Status = simwing::fsi::
            ScenePressureCellMimeticConductancePhaseSampleStatus;
        if (sample.status == Status::Accepted) {
            accepted_.push_back(sample.normalizedConductance);
            return;
        }
        if (sample.status == Status::RejectedLocalCellLinearConsistency) {
            ++rejectedCount_;
            return;
        }
        throw std::logic_error("mimetic conductance sample status is unknown");
    }

    void print() const {
        double minimum = 0.0;
        double maximum = 0.0;
        double mean = 0.0;
        double variation = 0.0;
        if (!accepted_.empty()) {
            minimum = *std::min_element(accepted_.begin(), accepted_.end());
            maximum = *std::max_element(accepted_.begin(), accepted_.end());
            CompensatedSum sum;
            for (const double value : accepted_) {
                sum.add(value);
            }
            mean = sum.value() / static_cast<double>(accepted_.size());
            CompensatedSum variance;
            for (const double value : accepted_) {
                const double difference = value - mean;
                variance.add(difference * difference);
            }
            variation = std::sqrt(std::max(
                0.0,
                variance.value() / static_cast<double>(accepted_.size())))
                / mean;
        }
        std::cout << "summary accepted=" << accepted_.size()
                  << " rejected=" << rejectedCount_
                  << " minimum=" << minimum
                  << " maximum=" << maximum
                  << " mean=" << mean
                  << " cv=" << variation << '\n';
    }

private:
    std::vector<double> accepted_;
    std::size_t rejectedCount_ = 0;
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
        "runs sequentially with one heavy immutable product live at a time. "
        "Fine-grid\n"
        "samples can take minutes. The tool\n"
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

const simwing::fsi::ScenePressureCellMimeticConductancePhaseSample&
printAuditSample(
    const simwing::fsi::
        ScenePressureCellMimeticConductancePhaseRefinementAudit& audit,
    const std::size_t canonicalPhaseIndex) {
    if (audit.levels.size() != 1
        || audit.levels.front().samples.size() != 1) {
        throw std::logic_error(
            "streamed mimetic conductance audit is not a single sample");
    }
    const auto& level = audit.levels.front();
    const auto& sample = level.samples.front();
    const auto phase = sample.gridPhaseFraction;
    std::cout << "sample phase=" << canonicalPhaseIndex
              << " translation=(" << phase.x << ',' << phase.y << ','
              << phase.z << ')'
              << " status=" << statusName(sample.status)
              << " audit-fingerprint=0x" << std::hex << audit.fingerprint
              << std::dec
              << " owned-bytes=" << audit.ownedStorageBytes
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
                  << rejection.localCell.maximumAlgebraicConsistencyError
                  << " algebraic-tolerance="
                  << rejection.localCell.algebraicConsistencyTolerance;
    }
    std::cout << std::endl;
    return sample;
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
        if (options.allPhases) {
            canonicalPhaseIndices.reserve(
                simwing::fsi::
                    scenePressureCellMimeticConductanceCanonicalGridPhases
                        .size());
            for (std::size_t index = 0;
                 index < simwing::fsi::
                     scenePressureCellMimeticConductanceCanonicalGridPhases
                         .size();
                 ++index) {
                canonicalPhaseIndices.push_back(index);
            }
        } else {
            canonicalPhaseIndices.push_back(*options.phaseIndex);
        }
        const std::array<simwing::fsi::fluid::GridCellCounts, 1>
            resolutions{{{
                options.resolution, options.resolution,
                options.resolution,
            }}};
        const auto settings = simwing::fsi::
            makeScenePressureCellMimeticConductanceOfflineAuditSettings();
        std::cout.imbue(std::locale::classic());
        std::cout << std::setprecision(
            std::numeric_limits<double>::max_digits10);
        std::cout << "run audit-version="
                  << simwing::fsi::
                         scenePressureCellMimeticConductancePhaseRefinementAuditVersion
                  << " resolution=" << options.resolution
                  << " phase-count=" << canonicalPhaseIndices.size()
                  << " execution=sequential" << std::endl;
        CompactAggregate aggregate;
        std::uint64_t structureFingerprint = 0;
        for (const std::size_t canonicalPhaseIndex
             : canonicalPhaseIndices) {
            const std::array<simwing::fsi::fluid::Vector3, 1> phases{{
                simwing::fsi::
                    scenePressureCellMimeticConductanceCanonicalGridPhases[
                        canonicalPhaseIndex],
            }};
            const auto audit = simwing::fsi::
                auditScenePressureCellMimeticConductancePhaseRefinement(
                    resolutions, phases, settings);
            if (structureFingerprint == 0) {
                structureFingerprint =
                    audit.structureDefinitionFingerprint;
            } else if (audit.structureDefinitionFingerprint
                       != structureFingerprint) {
                throw std::logic_error(
                    "streamed mimetic conductance structure changed");
            }
            aggregate.add(printAuditSample(audit, canonicalPhaseIndex));
        }
        std::cout << "structure-fingerprint=0x" << std::hex
                  << structureFingerprint << std::dec << '\n';
        aggregate.print();
        return 0;
    } catch (const std::exception& exception) {
        std::fprintf(stderr, "error: audit failed: %s\n", exception.what());
        return 1;
    }
}
