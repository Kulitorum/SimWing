foreach(required
        SIMWING_FSI
        CONTROL_FIXTURE
        COMMAND_FILE
        RESPONSE_FILE
        CHECKPOINT_FILE
        TRACE_FILE
        ERROR_FILE)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "${required} is required")
    endif()
endforeach()

file(REMOVE
    "${COMMAND_FILE}"
    "${RESPONSE_FILE}"
    "${CHECKPOINT_FILE}"
    "${TRACE_FILE}"
    "${ERROR_FILE}")

execute_process(
    COMMAND "${CONTROL_FIXTURE}" write "${COMMAND_FILE}"
    RESULT_VARIABLE write_result
    OUTPUT_VARIABLE write_output
    ERROR_VARIABLE write_error
    TIMEOUT 15
)
if(NOT write_result EQUAL 0)
    message(FATAL_ERROR
        "control command fixture failed (${write_result}):\n"
        "${write_output}${write_error}")
endif()

execute_process(
    COMMAND "${SIMWING_FSI}"
        --case periodic-flow
        --control-stdio
        --steps 1
        --no-viewer
    RESULT_VARIABLE invalid_steps_result
    OUTPUT_VARIABLE invalid_steps_output
    ERROR_VARIABLE invalid_steps_error
    TIMEOUT 15
)
if(invalid_steps_result EQUAL 0)
    message(FATAL_ERROR
        "control mode unexpectedly accepted --steps:\n"
        "${invalid_steps_output}${invalid_steps_error}")
endif()

execute_process(
    COMMAND "${SIMWING_FSI}"
        --case periodic-flow
        --control-stdio
        --trace "${TRACE_FILE}"
        --checkpoint-out "${CHECKPOINT_FILE}"
    INPUT_FILE "${COMMAND_FILE}"
    OUTPUT_FILE "${RESPONSE_FILE}"
    ERROR_FILE "${ERROR_FILE}"
    RESULT_VARIABLE worker_result
    TIMEOUT 30
)
if(NOT worker_result EQUAL 0)
    file(READ "${ERROR_FILE}" worker_error)
    message(FATAL_ERROR
        "controlled worker failed (${worker_result}):\n${worker_error}")
endif()

execute_process(
    COMMAND "${CONTROL_FIXTURE}" verify
        "${RESPONSE_FILE}"
        "${CHECKPOINT_FILE}"
        "${TRACE_FILE}"
    RESULT_VARIABLE verify_result
    OUTPUT_VARIABLE verify_output
    ERROR_VARIABLE verify_error
    TIMEOUT 15
)
if(NOT verify_result EQUAL 0)
    file(READ "${ERROR_FILE}" worker_error)
    message(FATAL_ERROR
        "control response fixture failed (${verify_result}):\n"
        "${verify_output}${verify_error}\nworker stderr:\n${worker_error}")
endif()
