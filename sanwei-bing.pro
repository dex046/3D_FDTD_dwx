TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    DataTran.cpp \
    H_Border.cpp \
    h_Coord.cpp \
    Inside.cpp \
    Partition.cpp \
    RWsgy.cpp \
    testHFWI3D.cpp

include(deployment.pri)
qtcAddDeployment()

INCLUDEPATH += /usr/include/mpich-x86_64/
DEPENDPATH  += /usr/include/mpich-x86_64/

INCLUDEPATH += /usr/local/fftw-3.3.4/include/
DEPENDPATH  += /usr/local/fftw-3.3.4/include/

HEADERS += \
    DataTran.h \
    Partition.h \
    RWsgy.h \
    testHFWI3D.h


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/release/ -lfmpich
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/debug/ -lfmpich
else:unix: LIBS += -L$$PWD/../../../lib64/mpich/lib/ -lfmpich

INCLUDEPATH += $$PWD/../../../lib64/mpich
DEPENDPATH += $$PWD/../../../lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/release/ -lmpich
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/debug/ -lmpich
else:unix: LIBS += -L$$PWD/../../../lib64/mpich/lib/ -lmpich

INCLUDEPATH += $$PWD/../../../lib64/mpich
DEPENDPATH += $$PWD/../../../lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/release/ -lmpichcxx
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/debug/ -lmpichcxx
else:unix: LIBS += -L$$PWD/../../../lib64/mpich/lib/ -lmpichcxx

INCLUDEPATH += $$PWD/../../../lib64/mpich
DEPENDPATH += $$PWD/../../../lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/release/ -lmpichf90
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/debug/ -lmpichf90
else:unix: LIBS += -L$$PWD/../../../lib64/mpich/lib/ -lmpichf90

INCLUDEPATH += $$PWD/../../../lib64/mpich
DEPENDPATH += $$PWD/../../../lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/release/ -lmpl
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/debug/ -lmpl
else:unix: LIBS += -L$$PWD/../../../lib64/mpich/lib/ -lmpl

INCLUDEPATH += $$PWD/../../../lib64/mpich
DEPENDPATH += $$PWD/../../../lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/release/ -lopa
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib64/mpich/lib/debug/ -lopa
else:unix: LIBS += -L$$PWD/../../../lib64/mpich/lib/ -lopa

INCLUDEPATH += $$PWD/../../../lib64/mpich
DEPENDPATH += $$PWD/../../../lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/release/ -lfftw3
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/debug/ -lfftw3
else:unix: LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/ -lfftw3

INCLUDEPATH += $$PWD/../../../usr/local/fftw-3.3.4/include
DEPENDPATH += $$PWD/../../../usr/local/fftw-3.3.4/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/release/ -lfmpich
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/debug/ -lfmpich
else:unix: LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/ -lfmpich

INCLUDEPATH += $$PWD/../../../usr/lib64/mpich
DEPENDPATH += $$PWD/../../../usr/lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/release/ -lmpich
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/debug/ -lmpich
else:unix: LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/ -lmpich

INCLUDEPATH += $$PWD/../../../usr/lib64/mpich
DEPENDPATH += $$PWD/../../../usr/lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/release/ -lmpichcxx
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/debug/ -lmpichcxx
else:unix: LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/ -lmpichcxx

INCLUDEPATH += $$PWD/../../../usr/lib64/mpich
DEPENDPATH += $$PWD/../../../usr/lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/release/ -lmpichf90
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/debug/ -lmpichf90
else:unix: LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/ -lmpichf90

INCLUDEPATH += $$PWD/../../../usr/lib64/mpich
DEPENDPATH += $$PWD/../../../usr/lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/release/ -lmpl
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/debug/ -lmpl
else:unix: LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/ -lmpl

INCLUDEPATH += $$PWD/../../../usr/lib64/mpich
DEPENDPATH += $$PWD/../../../usr/lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/release/ -lopa
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/debug/ -lopa
else:unix: LIBS += -L$$PWD/../../../usr/lib64/mpich/lib/ -lopa

INCLUDEPATH += $$PWD/../../../usr/lib64/mpich
DEPENDPATH += $$PWD/../../../usr/lib64/mpich

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/release/ -lfftw3
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/debug/ -lfftw3
else:unix: LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/ -lfftw3

INCLUDEPATH += $$PWD/../../../usr/local/fftw-3.3.4/include
DEPENDPATH += $$PWD/../../../usr/local/fftw-3.3.4/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/release/ -lfftw3
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/debug/ -lfftw3
else:unix: LIBS += -L$$PWD/../../../usr/local/fftw-3.3.4/lib/ -lfftw3

INCLUDEPATH += $$PWD/../../../usr/local/fftw-3.3.4/include
DEPENDPATH += $$PWD/../../../usr/local/fftw-3.3.4/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../usr/local/fftw-3.3.4/lib/release/libfftw3.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../usr/local/fftw-3.3.4/lib/debug/libfftw3.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../usr/local/fftw-3.3.4/lib/release/fftw3.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../usr/local/fftw-3.3.4/lib/debug/fftw3.lib
else:unix: PRE_TARGETDEPS += $$PWD/../../../usr/local/fftw-3.3.4/lib/libfftw3.a

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../usr/lib/fftw-3.3.4/release/ -lfftw3
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../usr/lib/fftw-3.3.4/debug/ -lfftw3
else:unix: LIBS += -L$$PWD/../../../usr/lib/fftw-3.3.4/ -lfftw3

INCLUDEPATH += $$PWD/../../../usr/lib/fftw-3.3.4
DEPENDPATH += $$PWD/../../../usr/lib/fftw-3.3.4

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../lib/fftw-3.3.4/release/ -lfftw3
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../lib/fftw-3.3.4/debug/ -lfftw3
else:unix: LIBS += -L$$PWD/../../../lib/fftw-3.3.4/ -lfftw3

INCLUDEPATH += $$PWD/../../../lib/fftw-3.3.4
DEPENDPATH += $$PWD/../../../lib/fftw-3.3.4

OTHER_FILES += \
    information.txt
