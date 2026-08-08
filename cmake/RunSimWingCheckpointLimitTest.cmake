if(NOT DEFINED SIMWING_FSI
   OR NOT DEFINED TRACE_FILE
   OR NOT DEFINED CHECKPOINT_FILE)
    message(FATAL_ERROR "checkpoint limit test inputs are incomplete")
endif()

file(REMOVE "${TRACE_FILE}" "${CHECKPOINT_FILE}")
execute_process(
    COMMAND "${SIMWING_FSI}"
        --case moving-porous-flow
        --steps 10001
        --no-viewer
        --trace "${TRACE_FILE}"
        --checkpoint-out "${CHECKPOINT_FILE}"
    RESULT_VARIABLE worker_result
    OUTPUT_VARIABLE worker_stdout
    ERROR_VARIABLE worker_stderr)

if(NOT worker_result EQUAL 2)
    message(FATAL_ERROR
        "worker returned ${worker_result}, expected option rejection 2\n"
        "stdout: ${worker_stdout}\n"
        "stderr: ${worker_stderr}")
endif()
if(NOT worker_stderr MATCHES
   "moving-porous-flow checkpoint output would exceed the 10000-step deterministic replay limit")
    message(FATAL_ERROR "unexpected worker diagnostic: ${worker_stderr}")
endif()
if(EXISTS "${TRACE_FILE}" OR EXISTS "${CHECKPOINT_FILE}")
    message(FATAL_ERROR
        "checkpoint replay-limit rejection created an output file")
endif()
