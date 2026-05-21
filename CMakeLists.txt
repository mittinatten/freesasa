cmake_minimum_required(VERSION 3.14)
project(FreeSASA VERSION 2.1.3 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Options ───────────────────────────────────────────────────────────────────
option(ENABLE_JSON   "Build with JSON output support (requires json-c)"   OFF)
option(ENABLE_XML    "Build with XML output support (requires libxml2)"   OFF)
option(ENABLE_TESTS  "Build the unit test suite (requires Check)"         ON)

# ── Math library ──────────────────────────────────────────────────────────────
find_library(MATH_LIBRARY m)

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
            "Install with: brew install check  OR  conda install -c conda-forge check")
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
    USE_THREADS=0             # threading off in base CMake build; add OpenMP separately
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
if(ENABLE_JSON)
    target_include_directories(freesasa_static PUBLIC ${JSON_C_INCLUDE_DIRS})
    target_link_libraries(freesasa_static PUBLIC ${JSON_C_LIBRARIES})
endif()
if(ENABLE_XML)
    target_include_directories(freesasa_static PUBLIC ${LIBXML2_INCLUDE_DIRS})
    target_link_libraries(freesasa_static PUBLIC ${LIBXML2_LIBRARIES})
endif()

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
if(ENABLE_JSON)
    target_include_directories(freesasa_shared PUBLIC ${JSON_C_INCLUDE_DIRS})
    target_link_libraries(freesasa_shared PUBLIC ${JSON_C_LIBRARIES})
endif()
if(ENABLE_XML)
    target_include_directories(freesasa_shared PUBLIC ${LIBXML2_INCLUDE_DIRS})
    target_link_libraries(freesasa_shared PUBLIC ${LIBXML2_LIBRARIES})
endif()

# ── CLI executable ────────────────────────────────────────────────────────────
add_executable(freesasa_bin src/main.cc src/cif.cc)
set_target_properties(freesasa_bin PROPERTIES OUTPUT_NAME freesasa)
target_compile_definitions(freesasa_bin PRIVATE ${FREESASA_COMPILE_DEFS})
target_include_directories(freesasa_bin PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)
target_link_libraries(freesasa_bin PRIVATE freesasa_static)

# ── Test suite ────────────────────────────────────────────────────────────────
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

    target_compile_definitions(freesasa_tests PRIVATE ${TEST_DEFS})
    target_include_directories(freesasa_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/tests
        ${CHECK_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/third-party/gemmi/include)

    find_package(Threads REQUIRED)
    target_link_libraries(freesasa_tests PRIVATE
        freesasa_test_static
        ${CHECK_LIBRARY}
        ${MATH_LIBRARY}
        Threads::Threads
        ${CMAKE_DL_LIBS})

    enable_testing()
    add_test(NAME freesasa_unit_tests
        COMMAND freesasa_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

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
message(STATUS "  Unit tests:   ${ENABLE_TESTS}")
message(STATUS "  C compiler:   ${CMAKE_C_COMPILER}")
message(STATUS "")
