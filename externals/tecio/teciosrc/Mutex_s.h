 #pragma once
#include "ThirdPartyHeadersBegin.h"
 #if defined _WIN32
#include <windows.h>
 #else
#include <pthread.h>
 #if defined __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
 #endif
 #endif
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
struct ___2665 {
 #if defined _WIN32
HANDLE ___2494;
 #else
pthread_mutex_t ___2494;
 #endif
___2665() {
 #if defined _WIN32
___2494 = CreateMutex(NULL, ___1305, NULL);
 #else
pthread_mutex_init(&___2494, NULL);
 #endif
} ~___2665() {
 #if defined _WIN32
CloseHandle(___2494);
 #else
pthread_mutex_destroy(&___2494);
 #endif
} void lock() {
 #if defined _WIN32
WaitForSingleObject(___2494, INFINITE);
 #else
pthread_mutex_lock(&___2494);
 #endif
} void unlock() {
 #if defined _WIN32
ReleaseMutex(___2494);
 #else
pthread_mutex_unlock(&___2494);
 #endif
} };
