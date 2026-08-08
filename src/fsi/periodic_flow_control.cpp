#include "periodic_flow_control.h"

#include <exception>
#include <utility>

namespace simwing::fsi {
namespace {

constexpr std::size_t maximumErrorMessageBytes =
    WorkerControlProtocolLimits{}.maximumErrorMessageBytes;

std::string boundedFailureMessage(
    std::string prefix,
    const std::string& detail) {
    if (!detail.empty()) {
        prefix += ": ";
        prefix += detail;
    }
    if (prefix.empty()) {
        prefix = "worker operation failed";
    }
    if (prefix.size() > maximumErrorMessageBytes) {
        prefix.resize(maximumErrorMessageBytes);
    }
    return prefix;
}

} // namespace

PeriodicFlowControlSession::PeriodicFlowControlSession(
    PeriodicFlowCase& simulation,
    PeriodicFlowControlHooks hooks)
    : simulation_(&simulation), hooks_(std::move(hooks)) {}

WorkerControlResponse PeriodicFlowControlSession::readyResponse() const {
    return response(WorkerControlResponseKind::Ready, 0);
}

bool PeriodicFlowControlSession::execute(
    const WorkerControlCommand& command,
    WorkerControlResponse& result,
    WorkerControlProtocolError* protocolError) {
    if (!validateWorkerControlCommand(command, protocolError)) {
        return false;
    }

    if (stopped_) {
        if (command.kind == WorkerControlCommandKind::Stop) {
            result = response(
                WorkerControlResponseKind::Stopped, command.requestId);
        } else {
            result = failureResponse(
                command.requestId,
                WorkerControlFailureCode::InvalidCommand,
                "worker has already stopped");
        }
        return true;
    }

    switch (command.kind) {
    case WorkerControlCommandKind::Advance: {
        std::uint64_t producedFrameCount = 0;
        for (std::uint64_t step = 0;
             step < command.advanceStepCount; ++step) {
            viewer::DiagnosticFrame frame;
            try {
                frame = simulation_->advance();
            } catch (const std::exception& exception) {
                result = failureResponse(
                    command.requestId,
                    WorkerControlFailureCode::NumericalFailure,
                    boundedFailureMessage(
                        "periodic-flow advance failed", exception.what()));
                return true;
            } catch (...) {
                result = failureResponse(
                    command.requestId,
                    WorkerControlFailureCode::NumericalFailure,
                    "periodic-flow advance failed with an unknown error");
                return true;
            }

            ++producedFrameCount;
            if (hooks_.publishFrame) {
                std::string hookError;
                try {
                    if (!hooks_.publishFrame(frame, hookError)) {
                        result = failureResponse(
                            command.requestId,
                            WorkerControlFailureCode::InternalFailure,
                            boundedFailureMessage(
                                "accepted frame publication failed",
                                hookError));
                        return true;
                    }
                } catch (const std::exception& exception) {
                    result = failureResponse(
                        command.requestId,
                        WorkerControlFailureCode::InternalFailure,
                        boundedFailureMessage(
                            "accepted frame publication threw",
                            exception.what()));
                    return true;
                } catch (...) {
                    result = failureResponse(
                        command.requestId,
                        WorkerControlFailureCode::InternalFailure,
                        "accepted frame publication threw an unknown error");
                    return true;
                }
            }
        }
        result = response(
            WorkerControlResponseKind::Advanced, command.requestId);
        result.producedFrameCount = producedFrameCount;
        return true;
    }
    case WorkerControlCommandKind::Checkpoint: {
        if (!hooks_.writeCheckpoint) {
            result = failureResponse(
                command.requestId,
                WorkerControlFailureCode::CheckpointFailure,
                "checkpoint persistence is not configured");
            return true;
        }
        try {
            const PeriodicFlowCaseCheckpoint checkpoint =
                simulation_->checkpoint();
            std::string hookError;
            if (!hooks_.writeCheckpoint(checkpoint, hookError)) {
                result = failureResponse(
                    command.requestId,
                    WorkerControlFailureCode::CheckpointFailure,
                    boundedFailureMessage(
                        "checkpoint persistence failed", hookError));
                return true;
            }
        } catch (const std::exception& exception) {
            result = failureResponse(
                command.requestId,
                WorkerControlFailureCode::CheckpointFailure,
                boundedFailureMessage(
                    "checkpoint persistence threw", exception.what()));
            return true;
        } catch (...) {
            result = failureResponse(
                command.requestId,
                WorkerControlFailureCode::CheckpointFailure,
                "checkpoint persistence threw an unknown error");
            return true;
        }
        result = response(
            WorkerControlResponseKind::Checkpointed, command.requestId);
        return true;
    }
    case WorkerControlCommandKind::Stop:
        stopped_ = true;
        result = response(
            WorkerControlResponseKind::Stopped, command.requestId);
        return true;
    }

    return false;
}

bool PeriodicFlowControlSession::stopped() const noexcept {
    return stopped_;
}

WorkerControlResponse PeriodicFlowControlSession::response(
    const WorkerControlResponseKind kind,
    const std::uint64_t requestId) const {
    WorkerControlResponse result;
    result.kind = kind;
    result.requestId = requestId;
    result.acceptedStepCount = simulation_->acceptedStepCount();
    result.simulationTimeSeconds = simulation_->simulationTimeSeconds();
    return result;
}

WorkerControlResponse PeriodicFlowControlSession::failureResponse(
    const std::uint64_t requestId,
    const WorkerControlFailureCode failureCode,
    std::string message) const {
    WorkerControlResponse result = response(
        WorkerControlResponseKind::Error, requestId);
    result.failureCode = failureCode;
    if (message.empty()) {
        message = "worker operation failed";
    }
    if (message.size() > maximumErrorMessageBytes) {
        message.resize(maximumErrorMessageBytes);
    }
    result.errorMessage = std::move(message);
    return result;
}

} // namespace simwing::fsi
