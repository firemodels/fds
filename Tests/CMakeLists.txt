
include(CTest)

set(FIG_COMMIT "513dc046131f5dfdda47ff485efc2db7243f06b1")

FetchContent_Declare(
    test_data
    GIT_REPOSITORY https://github.com/firemodels/fig.git
    GIT_TAG        ${FIG_COMMIT}
)
FetchContent_MakeAvailable(test_data)
set(FIG_DIR ${test_data_SOURCE_DIR})

add_test(NAME "Execute fds" COMMAND fds )
add_test(NAME "Execute fds JSON" COMMAND fds --json - ${FIG_DIR}/smv/Tests/test_model_small )
