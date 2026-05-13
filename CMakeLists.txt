cmake_minimum_required(VERSION 3.14)
project(FreeSASA VERSION 2.1.3 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ──────────────────────── Options ────────────────────────
option(ENABLE_JSON     "Build with JSON output support (requires json-c)"        OFF)
option(ENABLE_XML      "Build with XML output support (requires libxml2)"        OFF)
option(ENABLE_TESTS    "Build the Check-based test suite"                        ON)
option(ENABLE_MEMERR   "Include malloc-failure tests (requires dlsym hack)"      OFF)
option(ENABLE_PYTHON   "Install shared library for Python bindings"              ON)

# ──────────────────────── OpenMP ────────────────────────
find_package(OpenMP REQUIRED)
if(OpenMP_C_FOUND)
    message(STATUS "OpenMP C flags: ${OpenMP_C_FLAGS}")
endif()

# ──────────────────────── Math library ────────────────────────
find_library(MATH_LIBRARY m)

# ──────────────────────── Check (unit tests) ────────────────────────
# Look in conda/miniconda first, then system paths
set(CHECK_SEARCH_PATHS
    $ENV{HOME}/miniconda3
    $ENV{HOME}/miniconda
    $ENV{HOME}/anaconda3
    /usr/local
    /usr
)

find_path(CHECK_INCLUDE_DIR check.h
    PATHS ${CHECK_SEARCH_PATHS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH)
find_path(CHECK_INCLUDE_DIR check.h)   # fallback to default search

find_library(CHECK_LIBRARY
    NAMES check
    PATHS ${CHECK_SEARCH_PATHS}
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH)
find_library(CHECK_LIBRARY NAMES check) # fallback

if(CHECK_INCLUDE_DIR AND CHECK_LIBRARY)
    set(CHECK_FOUND TRUE)
    message(STATUS "Found Check: ${CHECK_LIBRARY}")
    message(STATUS "Check include: ${CHECK_INCLUDE_DIR}")
else()
    set(CHECK_FOUND FALSE)
    if(ENABLE_TESTS)
        message(WARNING "Check unit test framework not found — tests disabled. Install with: conda install -c conda-forge check")
        set(ENABLE_TESTS OFF)
    endif()
endif()

# ──────────────────────── JSON ────────────────────────
if(ENABLE_JSON)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSON_C REQUIRED json-c)
    set(USE_JSON 1)
else()
    set(USE_JSON 0)
endif()

# ──────────────────────── XML ────────────────────────
if(ENABLE_XML)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(LIBXML2 REQUIRED libxml-2.0)
    set(USE_XML 1)
else()
    set(USE_XML 0)
endif()

# ──────────────────────── Package defines (global) ────────────────────────
add_definitions(
    -DPACKAGE_VERSION="${PROJECT_VERSION}"
    -DPACKAGE_STRING="FreeSASA ${PROJECT_VERSION}"
    -DHAVE_CONFIG_H=1
    -DUSE_OPENMP=1
    -DUSE_THREADS=1
    -DUSE_JSON=${USE_JSON}
    -DUSE_XML=${USE_XML}
    -DUSE_CHECK=0
    -DREPORTBUG="Report bugs to <https://github.com/mittinatten/freesasa/issues>"
    -DHOMEPAGE="<http://freesasa.github.io>"
    -DFREESASA_XMLNS="http://freesasa.github.io/"
    -DSTDC_HEADERS=1
)

# ──────────────────────── Algorithm defaults ────────────────────────
add_definitions(
    -DFREESASA_DEF_ALGORITHM=FREESASA_LEE_RICHARDS
    -DFREESASA_DEF_PROBE_RADIUS=1.4
    -DFREESASA_DEF_SR_N=100
    -DFREESASA_DEF_LR_N=20
)

# ──────────────────────── Include paths ────────────────────────
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/third-party/gemmi/include
)

# ──────────────────────── Library sources ────────────────────────
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

# ──────────────────────── Static library ────────────────────────
add_library(freesasa_static STATIC ${LIBFREESASA_SOURCES})
set_target_properties(freesasa_static PROPERTIES OUTPUT_NAME freesasa)
target_link_libraries(freesasa_static PUBLIC OpenMP::OpenMP_C)
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

# Per-file SIMD flags for sasa_sr.c — enable AVX-512/AVX2 + fast-math
# These apply to all targets that compile sasa_sr.c
include(CheckCCompilerFlag)
check_c_compiler_flag("-march=native" HAS_MARCH_NATIVE)
check_c_compiler_flag("-ffast-math" HAS_FAST_MATH)
check_c_compiler_flag("-mavx512f" HAS_AVX512F)
check_c_compiler_flag("-mavx512dq" HAS_AVX512DQ)
check_c_compiler_flag("-mavx2" HAS_AVX2)
check_c_compiler_flag("-mfma" HAS_FMA)

set(SASA_SR_FLAGS "")
if(HAS_MARCH_NATIVE)
    list(APPEND SASA_SR_FLAGS "-march=native")
else()
    if(HAS_AVX512F AND HAS_AVX512DQ)
        list(APPEND SASA_SR_FLAGS "-mavx512f" "-mavx512dq")
    elseif(HAS_AVX2)
        list(APPEND SASA_SR_FLAGS "-mavx2")
    endif()
    if(HAS_FMA)
        list(APPEND SASA_SR_FLAGS "-mfma")
    endif()
endif()
if(HAS_FAST_MATH)
    list(APPEND SASA_SR_FLAGS "-ffast-math")
endif()
list(APPEND SASA_SR_FLAGS "-O3")

if(SASA_SR_FLAGS)
    set_source_files_properties(src/sasa_sr.c PROPERTIES
        COMPILE_OPTIONS "${SASA_SR_FLAGS}"
    )
    message(STATUS "sasa_sr.c SIMD flags: ${SASA_SR_FLAGS}")
endif()

# ──────────────────────── Shared library ────────────────────────
add_library(freesasa_shared SHARED ${LIBFREESASA_SOURCES})
set_target_properties(freesasa_shared PROPERTIES
    OUTPUT_NAME freesasa
    VERSION ${PROJECT_VERSION}
    SOVERSION 2
)
target_link_libraries(freesasa_shared PUBLIC OpenMP::OpenMP_C)
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

# ──────────────────────── CLI executable ────────────────────────
add_executable(freesasa_bin src/main.cc src/cif.cc)
set_target_properties(freesasa_bin PROPERTIES OUTPUT_NAME freesasa)
target_link_libraries(freesasa_bin PRIVATE freesasa_static OpenMP::OpenMP_CXX)

# ──────────────────────── Test suite ────────────────────────
if(ENABLE_TESTS AND CHECK_FOUND)
    # Absolute paths for test data and share directories
    set(TEST_DATADIR "${CMAKE_SOURCE_DIR}/tests/data/")
    set(TEST_SHAREDIR "${CMAKE_SOURCE_DIR}/share/")

    # Five source files have #if USE_CHECK embedded test functions:
    # classifier.c, nb.c, pdb.c, sasa_lr.c, selection.c
    # Recompile them with USE_CHECK=1 for the test binary only.
    set(CHECK_STUB_SOURCES
        src/classifier.c
        src/nb.c
        src/pdb.c
        src/sasa_lr.c
        src/selection.c
    )
    set(TEST_COMPILE_DEFS
        USE_CHECK=1
        USE_OPENMP=1
        USE_THREADS=1
        USE_JSON=${USE_JSON}
        USE_XML=${USE_XML}
        PACKAGE_VERSION="${PROJECT_VERSION}"
        PACKAGE_STRING="FreeSASA ${PROJECT_VERSION}"
        HAVE_CONFIG_H=1
        STDC_HEADERS=1
        DATADIR="${TEST_DATADIR}"
        SHAREDIR="${TEST_SHAREDIR}"
        INCLUDE_MEMERR_TESTS=$<IF:$<BOOL:${ENABLE_MEMERR}>,1,0>
        FREESASA_DEF_ALGORITHM=FREESASA_LEE_RICHARDS
        FREESASA_DEF_PROBE_RADIUS=1.4
        FREESASA_DEF_SR_N=100
        FREESASA_DEF_LR_N=20
        REPORTBUG="Report bugs to <https://github.com/mittinatten/freesasa/issues>"
        HOMEPAGE="<http://freesasa.github.io>"
        FREESASA_XMLNS="http://freesasa.github.io/"
    )

    add_library(freesasa_check_objs OBJECT ${CHECK_STUB_SOURCES})
    target_compile_definitions(freesasa_check_objs PRIVATE ${TEST_COMPILE_DEFS})
    target_include_directories(freesasa_check_objs PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}
        ${CHECK_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/third-party/gemmi/include
    )
    target_link_libraries(freesasa_check_objs PRIVATE OpenMP::OpenMP_C)

    set(TEST_SOURCES
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
    )

    add_executable(freesasa_tests
        ${TEST_SOURCES}
        $<TARGET_OBJECTS:freesasa_check_objs>
    )

    target_compile_definitions(freesasa_tests PRIVATE ${TEST_COMPILE_DEFS})

    target_include_directories(freesasa_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/tests
        ${CHECK_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/third-party/gemmi/include
    )

    # Build a reduced static lib that excludes the 5 files recompiled above
    set(LIBFREESASA_TEST_SOURCES ${LIBFREESASA_SOURCES})
    list(REMOVE_ITEM LIBFREESASA_TEST_SOURCES
        src/classifier.c
        src/nb.c
        src/pdb.c
        src/sasa_lr.c
        src/selection.c
    )

    add_library(freesasa_test_static STATIC ${LIBFREESASA_TEST_SOURCES})
    target_compile_definitions(freesasa_test_static PRIVATE
        USE_CHECK=0
        USE_OPENMP=1
        USE_THREADS=1
        USE_JSON=${USE_JSON}
        USE_XML=${USE_XML}
        PACKAGE_VERSION="${PROJECT_VERSION}"
        PACKAGE_STRING="FreeSASA ${PROJECT_VERSION}"
        HAVE_CONFIG_H=1
        STDC_HEADERS=1
        FREESASA_DEF_ALGORITHM=FREESASA_LEE_RICHARDS
        FREESASA_DEF_PROBE_RADIUS=1.4
        FREESASA_DEF_SR_N=100
        FREESASA_DEF_LR_N=20
        REPORTBUG="Report bugs to <https://github.com/mittinatten/freesasa/issues>"
        HOMEPAGE="<http://freesasa.github.io>"
        FREESASA_XMLNS="http://freesasa.github.io/"
    )
    target_include_directories(freesasa_test_static PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/third-party/gemmi/include
    )
    target_link_libraries(freesasa_test_static PUBLIC OpenMP::OpenMP_C)
    if(MATH_LIBRARY)
        target_link_libraries(freesasa_test_static PUBLIC ${MATH_LIBRARY})
    endif()

    target_link_libraries(freesasa_tests PRIVATE
        freesasa_test_static
        OpenMP::OpenMP_C
        ${CHECK_LIBRARY}
        ${MATH_LIBRARY}
    )

    # dlsym for malloc interception in tools.c
    target_link_libraries(freesasa_tests PRIVATE ${CMAKE_DL_LIBS})

    # pthread needed by Check
    find_package(Threads REQUIRED)
    target_link_libraries(freesasa_tests PRIVATE Threads::Threads)

    enable_testing()
    add_test(NAME freesasa_unit_tests
        COMMAND freesasa_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    message(STATUS "Tests: ENABLED (Check ${CHECK_LIBRARY})")
else()
    message(STATUS "Tests: DISABLED")
endif()

# ──────────────────────── Install ────────────────────────
include(GNUInstallDirs)

install(TARGETS freesasa_bin freesasa_static freesasa_shared
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES src/freesasa.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install classification data (used by Python bindings and CLI)
install(DIRECTORY share/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/freesasa
    FILES_MATCHING PATTERN "*.config"
)

# ──────────────────────── Summary ────────────────────────
message(STATUS "")
message(STATUS "FreeSASA ${PROJECT_VERSION} — Parallelized Build Configuration")
message(STATUS "  OpenMP:     ENABLED (runtime thread detection)")
message(STATUS "  JSON:       ${ENABLE_JSON}")
message(STATUS "  XML:        ${ENABLE_XML}")
message(STATUS "  Tests:      ${ENABLE_TESTS}")
message(STATUS "  MemErr:     ${ENABLE_MEMERR}")
message(STATUS "  C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "")
