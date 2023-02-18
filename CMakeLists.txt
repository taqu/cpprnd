cmake_minimum_required(VERSION 3.23)
project(cpprnd)

set(CMAKE_CXX_STANDARD 17)
set(SOURCES main.cpp cpprnd.h cpprnd.cpp)
add_executable(cpprnd ${SOURCES})

if(MSVC)

    set(DEFAULT_CXX_FLAGS "/DWIN32 /D_WINDOWS /D_UNICODE /DUNICODE /W4 /WX- /nologo /fp:precise /arch:AVX /Oi /Zc:wchar_t /TP /Gd")

    if(DEFINED USE_D3D12DOWNLEVEL)
        set(DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} /DLG3_D3D12_DOWNLEVEL")
    endif()

    if(MSVC_VERSION VERSION_LESS_EQUAL "1900")
        set(DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} /Zc:__cplusplus /std:c++latest")
    else()
        set(DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} /Zc:__cplusplus /std:c++17")
    endif()

    set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1 /Gy /GR- /GS /Gm- /EHsc")
    set(CMAKE_CXX_FLAGS_RELEASE "/MT /O2 /Oi /GL /GR- /DNDEBUG /EHsc-")

    set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "/LTCG")
elseif(UNIX)
    set(DEFAULT_CXX_FLAGS "-Wall -Wextra -O0 -std=c++17 -march=x86-64-v3 -fno-exceptions")
elseif(XCODE)
endif()
