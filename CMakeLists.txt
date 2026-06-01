cmake_minimum_required(VERSION 3.14)
project(FreeSASA VERSION 2.1.3 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)   # codebase uses C99 features (restrict, // comments)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Options ───────────────────────────────────────────────────────────────────
option(ENABLE_JSON              "Build with JSON output support (requires json-c)"   OFF)
option(ENABLE_XML               "Build with XML output support (requires libxml2)"   OFF)
option(ENABLE_THREADS           "Build with pthreads support"                        ON)
option(ENABLE_TESTS             "Build the unit test suite (requires Check)"         ON)
option(ENABLE_DOXYGEN           "Build API documentation with Doxygen"               OFF)
option(ENABLE_PARSER_GENERATOR  "Regenerate lexer/parser with Flex and Bison"        OFF)
option(ENABLE_GCOV              "Enable coverage instrumentation (gcov)"             OFF)

# ── Math library ──────────────────────────────────────────────────────────────
find_library(MATH_LIBRARY m)

# ── Threads (pthreads) ────────────────────────────────────────────────────────
if(ENABLE_THREADS)
    find_package(Threads REQUIRED)
    set(USE_THREADS 1)
else()
    set(USE_THREADS 0)
endif()

# ── Optional: JSON ─────────────────────────────────────────────────────────────
if(ENABLE_JSON)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSON_C REQUIRED json-c)
    set(USE_JSON 1)
else()
    set(USE_JSON 0)
endif()

# ── Optional: XML ──────────────────────────────────────────────────────────────
if(ENABLE_XML)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(LIBXML2 REQUIRED libxml-2.0)
    set(USE_XML 1)
else()
    set(USE_XML 0)
endif()

# ── Optional: Doxygen ─────────────────────────────────────────────────────────
if(ENABLE_DOXYGEN)
    find_package(Doxygen REQUIRED)
    set(PACKAGE_NAME    "${PROJECT_NAME}")
    set(PACKAGE_VERSION "${PROJECT_VERSION}")
    set(top_srcdir      "${CMAKE_SOURCE_DIR}")
    configure_file(
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${CMAKE_BINARY_DIR}/doc/Doxyfile
        @ONLY)
    # Doxyfile.in uses relative paths (doxygen-header.html, html/doxygen, etc.)
    # so doxygen must run from the doc/ directory for them to resolve correctly.
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
        COMMENT "Generating API documentation (cmake --build <dir> -t doc)"
        VERBATIM)
endif()

# ── Optional: Parser generator (Flex/Bison) ───────────────────────────────────
if(ENABLE_PARSER_GENERATOR)
    # Bison 3+ required; macOS ships Bison 2.3 which is too old.
    # Use -DBISON_ROOT=/opt/homebrew/opt/bison (or wherever Bison 3 is)
    # to point CMake at the right binary.
    find_package(FLEX REQUIRED)
    if(BISON_ROOT)
        find_package(BISON 3.0 REQUIRED PATHS ${BISON_ROOT} NO_DEFAULT_PATH)
    else()
        find_package(BISON 3.0 REQUIRED)
    endif()
    flex_target(FreesasaLexer
        ${CMAKE_SOURCE_DIR}/src/lexer.l
        ${CMAKE_SOURCE_DIR}/src/lexer.c)
    bison_target(FreesasaParser
        ${CMAKE_SOURCE_DIR}/src/parser.y
        ${CMAKE_SOURCE_DIR}/src/parser.c
        COMPILE_FLAGS "-d")
    add_flex_bison_dependency(FreesasaLexer FreesasaParser)
endif()

# ── Optional: gcov coverage ──────────────────────────────────────────────────
if(ENABLE_GCOV)
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   --coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

# ── Optional: Check (unit tests) ──────────────────────────────────────────────
# Search conda/miniconda/homebrew prefixes before falling back to system paths
set(FREESASA_HINT_PATHS
    "$ENV{HOME}/miniconda3"
    "$ENV{HOME}/miniconda"
    "$ENV{HOME}/anaconda3"
    "$ENV{HOME}/opt/homebrew"   # Apple Silicon Homebrew
    "/opt/homebrew"             # Apple Silicon Homebrew
    "/usr/local"                # Intel macOS Homebrew / Linux
    "/usr"
)

find_path(CHECK_INCLUDE_DIR check.h
    HINTS ${FREESASA_HINT_PATHS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH)
find_path(CHECK_INCLUDE_DIR check.h)   # system fallback

find_library(CHECK_LIBRARY
    NAMES check
    HINTS ${FREESASA_HINT_PATHS}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH)
find_library(CHECK_LIBRARY NAMES check)  # system fallback

if(CHECK_INCLUDE_DIR AND CHECK_LIBRARY)
    set(CHECK_FOUND TRUE)
    message(STATUS "Found Check: ${CHECK_LIBRARY}")
else()
    set(CHECK_FOUND FALSE)
    if(ENABLE_TESTS)
        message(WARNING
            "Check not found — tests disabled. "
            "Install with:  apt install check  |  brew install check  |  conda install -c conda-forge check")
        set(ENABLE_TESTS OFF)
    endif()
endif()

# ── Generate config.h in the source root (mirrors autotools ./configure output)
# Tests use #include "../config.h" with a relative path, so it must live here.
set(TEST_DATADIR  "${CMAKE_SOURCE_DIR}/tests/data/")
set(TEST_SHAREDIR "${CMAKE_SOURCE_DIR}/share/")
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/config.h.cmake
    ${CMAKE_SOURCE_DIR}/config.h
    @ONLY)
message(STATUS "Generated config.h in source root")

# ── Compile definitions shared by all targets ─────────────────────────────────
set(FREESASA_COMPILE_DEFS
    HAVE_CONFIG_H=1
    USE_THREADS=${USE_THREADS}
    USE_JSON=${USE_JSON}
    USE_XML=${USE_XML}
    FREESASA_DEF_ALGORITHM=FREESASA_LEE_RICHARDS
    FREESASA_DEF_PROBE_RADIUS=1.4
    FREESASA_DEF_SR_N=100
    FREESASA_DEF_LR_N=20
)

# ── Library sources ───────────────────────────────────────────────────────────
set(LIBFREESASA_SOURCES
    src/classifier.c
    src/classifier_protor.c
    src/classifier_oons.c
    src/classifier_naccess.c
    src/coord.c
    src/pdb.c
    src/log.c
    src/sasa_lr.c
    src/sasa_sr.c
    src/structure.c
    src/node.c
    src/freesasa.c
    src/nb.c
    src/util.c
    src/rsa.c
    src/selection.c
    src/lexer.c
    src/parser.c
)

if(ENABLE_JSON)
    list(APPEND LIBFREESASA_SOURCES src/json.c)
endif()
if(ENABLE_XML)
    list(APPEND LIBFREESASA_SOURCES src/xml.c)
endif()

# ── Helper: apply JSON/XML link flags to a target ─────────────────────────────
# pkg_check_modules sets _INCLUDE_DIRS, _LIBRARY_DIRS and _LINK_LIBRARIES.
# We need _LIBRARY_DIRS for the linker to find the .dylib/.so on macOS/Linux.
function(freesasa_link_optional_deps target)
    if(ENABLE_JSON)
        target_include_directories(${target} PUBLIC ${JSON_C_INCLUDE_DIRS})
        target_link_directories(${target}    PUBLIC ${JSON_C_LIBRARY_DIRS})
        target_link_libraries(${target}      PUBLIC ${JSON_C_LINK_LIBRARIES})
    endif()
    if(ENABLE_XML)
        target_include_directories(${target} PUBLIC ${LIBXML2_INCLUDE_DIRS})
        target_link_directories(${target}    PUBLIC ${LIBXML2_LIBRARY_DIRS})
        target_link_libraries(${target}      PUBLIC ${LIBXML2_LINK_LIBRARIES})
    endif()
    if(ENABLE_THREADS)
        target_link_libraries(${target} PUBLIC Threads::Threads)
    endif()
endfunction()

# ── Static library ────────────────────────────────────────────────────────────
add_library(freesasa_static STATIC ${LIBFREESASA_SOURCES})
set_target_properties(freesasa_static PROPERTIES OUTPUT_NAME freesasa)
target_compile_definitions(freesasa_static PRIVATE ${FREESASA_COMPILE_DEFS})
target_include_directories(freesasa_static PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)
if(MATH_LIBRARY)
    target_link_libraries(freesasa_static PUBLIC ${MATH_LIBRARY})
endif()
freesasa_link_optional_deps(freesasa_static)

# ── Shared library ────────────────────────────────────────────────────────────
add_library(freesasa_shared SHARED ${LIBFREESASA_SOURCES})
set_target_properties(freesasa_shared PROPERTIES
    OUTPUT_NAME freesasa
    VERSION     ${PROJECT_VERSION}
    SOVERSION   2)
target_compile_definitions(freesasa_shared PRIVATE ${FREESASA_COMPILE_DEFS})
target_include_directories(freesasa_shared PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)
if(MATH_LIBRARY)
    target_link_libraries(freesasa_shared PUBLIC ${MATH_LIBRARY})
endif()
freesasa_link_optional_deps(freesasa_shared)

# ── CLI executable ────────────────────────────────────────────────────────────
add_executable(freesasa_bin src/main.cc src/cif.cc)
set_target_properties(freesasa_bin PROPERTIES OUTPUT_NAME freesasa)
target_compile_definitions(freesasa_bin PRIVATE ${FREESASA_COMPILE_DEFS})
target_include_directories(freesasa_bin PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)
target_link_libraries(freesasa_bin PRIVATE freesasa_static)

# ── Test suite ────────────────────────────────────────────────────────────────
enable_testing()

if(ENABLE_TESTS AND CHECK_FOUND)
    # TEST_DATADIR and TEST_SHAREDIR are set above for configure_file()

    # Five source files embed test functions behind #if USE_CHECK.
    # Recompile them with USE_CHECK=1 for the test binary only.
    set(CHECK_STUB_SOURCES
        src/classifier.c
        src/nb.c
        src/pdb.c
        src/sasa_lr.c
        src/selection.c)

    set(TEST_DEFS
        ${FREESASA_COMPILE_DEFS}
        USE_CHECK=1
        DATADIR="${TEST_DATADIR}"
        SHAREDIR="${TEST_SHAREDIR}"
        INCLUDE_MEMERR_TESTS=0)

    add_library(freesasa_check_objs OBJECT ${CHECK_STUB_SOURCES})
    target_compile_definitions(freesasa_check_objs PRIVATE ${TEST_DEFS})
    target_include_directories(freesasa_check_objs PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}
        ${CHECK_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)

    # Reduced static lib excluding the five files recompiled above
    set(LIBFREESASA_TEST_SOURCES ${LIBFREESASA_SOURCES})
    list(REMOVE_ITEM LIBFREESASA_TEST_SOURCES
        src/classifier.c src/nb.c src/pdb.c src/sasa_lr.c src/selection.c)

    add_library(freesasa_test_static STATIC ${LIBFREESASA_TEST_SOURCES})
    target_compile_definitions(freesasa_test_static PRIVATE ${FREESASA_COMPILE_DEFS})
    target_include_directories(freesasa_test_static PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)
    if(MATH_LIBRARY)
        target_link_libraries(freesasa_test_static PUBLIC ${MATH_LIBRARY})
    endif()
    freesasa_link_optional_deps(freesasa_test_static)

    add_executable(freesasa_tests
        tests/test_main.c
        tests/test_freesasa.c
        tests/test_nb.c
        tests/test_pdb.c
        tests/test_coord.c
        tests/test_classifier.c
        tests/test_structure.c
        tests/test_node.c
        tests/test_selection.c
        tests/tools.c
        $<TARGET_OBJECTS:freesasa_check_objs>)

    # test_json.c and test_xml.c are only compiled when those features are on,
    # matching the autotools conditional (test_api_SOURCES += test_json.c etc.)
    if(ENABLE_JSON)
        target_sources(freesasa_tests PRIVATE tests/test_json.c)
    endif()
    if(ENABLE_XML)
        target_sources(freesasa_tests PRIVATE tests/test_xml.c)
    endif()

    target_compile_definitions(freesasa_tests PRIVATE ${TEST_DEFS})
    target_include_directories(freesasa_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/tests
        ${CHECK_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)

    find_package(Threads REQUIRED)

    # Ubuntu's libcheck is built with subunit support and requires -lsubunit.
    # Other platforms (macOS homebrew, conda) do not ship it — find optionally.
    find_library(SUBUNIT_LIBRARY NAMES subunit)

    target_link_libraries(freesasa_tests PRIVATE
        freesasa_test_static
        ${CHECK_LIBRARY}
        ${MATH_LIBRARY}
        Threads::Threads
        ${CMAKE_DL_LIBS})

    if(SUBUNIT_LIBRARY)
        target_link_libraries(freesasa_tests PRIVATE ${SUBUNIT_LIBRARY})
    endif()

    # When JSON/XML are enabled the test binary also links those libs directly
    # (same as autotools: test_api_LDADD += -ljson-c / ${libxml2_LIBS})
    if(ENABLE_JSON)
        target_link_directories(freesasa_tests PRIVATE ${JSON_C_LIBRARY_DIRS})
        target_link_libraries(freesasa_tests PRIVATE ${JSON_C_LINK_LIBRARIES})
    endif()
    if(ENABLE_XML)
        target_link_directories(freesasa_tests PRIVATE ${LIBXML2_LIBRARY_DIRS})
        target_link_libraries(freesasa_tests PRIVATE ${LIBXML2_LINK_LIBRARIES})
    endif()

    # Unit tests (test-api)
    add_test(NAME test-api
        COMMAND freesasa_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    # ── CLI test suite (test-cli) ─────────────────────────────────────────
    # Substitute autotools-style @VAR@ placeholders in test-cli.in
    set(top_srcdir "${CMAKE_SOURCE_DIR}")
    if(ENABLE_XML)
        set(USE_XML_STR "yes")
    else()
        set(USE_XML_STR "no")
    endif()
    if(ENABLE_JSON)
        set(USE_JSON_STR "yes")
    else()
        set(USE_JSON_STR "no")
    endif()
    find_program(XMLLINT_PROG xmllint)
    if(XMLLINT_PROG)
        set(XMLLINT "xmllint")
    else()
        set(XMLLINT "")
    endif()
    find_program(JSONLINT_PROG jsonlint)
    if(JSONLINT_PROG)
        set(JSONLINT "jsonlint")
    else()
        set(JSONLINT "")
    endif()

    # test-cli.in uses @USE_XML@ and @USE_JSON@ (not _STR), so alias them
    set(USE_XML_CONF  "${USE_XML_STR}")
    set(USE_JSON_CONF "${USE_JSON_STR}")

    # Read template file, substitute, write output
    # We do manual substitution because the template uses @USE_XML@ not @USE_XML_STR@
    file(READ ${CMAKE_SOURCE_DIR}/tests/test-cli.in TEST_CLI_CONTENT)
    string(REPLACE "@top_srcdir@" "${CMAKE_SOURCE_DIR}"  TEST_CLI_CONTENT "${TEST_CLI_CONTENT}")
    string(REPLACE "@USE_XML@"    "${USE_XML_STR}"       TEST_CLI_CONTENT "${TEST_CLI_CONTENT}")
    string(REPLACE "@USE_JSON@"   "${USE_JSON_STR}"      TEST_CLI_CONTENT "${TEST_CLI_CONTENT}")
    string(REPLACE "@XMLLINT@"    "${XMLLINT}"           TEST_CLI_CONTENT "${TEST_CLI_CONTENT}")
    string(REPLACE "@JSONLINT@"   "${JSONLINT}"          TEST_CLI_CONTENT "${TEST_CLI_CONTENT}")
    file(WRITE ${CMAKE_BINARY_DIR}/tests/test-cli "${TEST_CLI_CONTENT}")
    file(CHMOD ${CMAKE_BINARY_DIR}/tests/test-cli
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                     GROUP_READ GROUP_EXECUTE
                     WORLD_READ WORLD_EXECUTE)

    # Create tmp/ directory needed by the test script
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/tmp)

    # The test-cli script expects the freesasa binary at ../src/freesasa
    # relative to the tests/ working directory.  Create a symlink so
    # the CMake build layout matches what the script expects.
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/src)
    add_custom_command(TARGET freesasa_bin POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE:freesasa_bin>
            ${CMAKE_BINARY_DIR}/src/freesasa
        COMMENT "Symlinking freesasa binary for CLI tests")

    add_test(NAME test-cli
        COMMAND bash ${CMAKE_BINARY_DIR}/tests/test-cli
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

    message(STATUS "Tests: ENABLED (Check: ${CHECK_LIBRARY})")
else()
    message(STATUS "Tests: DISABLED")
endif()

# ── Install ───────────────────────────────────────────────────────────────────
include(GNUInstallDirs)
install(TARGETS freesasa_bin freesasa_static freesasa_shared
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES src/freesasa.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY share/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/freesasa
    FILES_MATCHING PATTERN "*.config")

# ── Summary ───────────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "FreeSASA ${PROJECT_VERSION} — CMake configuration summary")
message(STATUS "  JSON output:  ${ENABLE_JSON}")
message(STATUS "  XML output:   ${ENABLE_XML}")
message(STATUS "  Threads:      ${ENABLE_THREADS}")
message(STATUS "  Unit tests:   ${ENABLE_TESTS}")
message(STATUS "  Doxygen:      ${ENABLE_DOXYGEN}")
message(STATUS "  Parser gen:   ${ENABLE_PARSER_GENERATOR}")
message(STATUS "  Coverage:     ${ENABLE_GCOV}")
message(STATUS "  C compiler:   ${CMAKE_C_COMPILER}")
message(STATUS "")
