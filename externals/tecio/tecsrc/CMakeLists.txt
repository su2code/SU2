cmake_minimum_required (VERSION 2.8)

project (tecio)

add_library(
    tecio SHARED
    alloc.cpp
    ALLOC.h
    arrlist.cpp
    ARRLIST.h
    auxdata.cpp
    AUXDATA.h
    dataio.cpp
    DATAIO.h
    dataio4.cpp
    DATAIO4.h
    dataset.cpp
    DATASET.h
    dataset0.cpp
    DATASET0.h
    DATASHR.h
    datautil.cpp
    DATAUTIL.h
    FACE.h
    filestream.cpp
    FILESTREAM.h
    GEOM.h
    geom2.cpp
    GEOM2.h
    GLOBAL.h
    INPUT.h
    MASTER.h
    q_msg.cpp
    Q_MSG.h
    q_unicode.cpp
    Q_UNICODE.h
    set.cpp
    SET.h
    stdafx.h
    strlist.cpp
    STRLIST.h
    strutil.cpp
    STRUTIL.h
    SYSTEM.h
    tassert.cpp
    TASSERT.h
    TECIO.h
    tecxxx.cpp
    TECXXX.h
    TEXT.h
    TranslatedString.cpp
    TranslatedString.h
    VERSION.h
    W__BASE.h
    )

set_target_properties(
    tecio PROPERTIES
    COMPILE_DEFINITIONS MAKEARCHIVE
    COMPILE_DEFINITIONS_RELEASE NO_ASSERTS
    )