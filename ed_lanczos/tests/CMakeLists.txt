# No automated tests exist yet. This file is a placeholder so that
# `-DED_LANCZOS_BUILD_TESTS=ON` doesn't fail the configure step.
#
# When tests are added (e.g. with Catch2 or GoogleTest), register them here
# with add_executable()/add_test(), for example:
#
#   find_package(Catch2 REQUIRED)
#   add_executable(unit_tests test_common.cpp test_kspace_1o2.cpp)
#   target_link_libraries(unit_tests PRIVATE ed_lanczos::ed_lanczos Catch2::Catch2WithMain)
#   add_test(NAME unit_tests COMMAND unit_tests)
message(STATUS "ed_lanczos: no unit tests defined yet")
