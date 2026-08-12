
# Fetch Google Test from GitHub
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest
  GIT_TAG        v1.15.2
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Add your test executable
add_executable(my_tests test_main.cpp)

# Link against Google Test (gtest_main provides a default main() function)
target_link_libraries(my_tests PRIVATE ed_lanczos::ed_lanczos PRIVATE gtest_main)

# Discover and run tests
include(GoogleTest)
gtest_discover_tests(my_tests)

