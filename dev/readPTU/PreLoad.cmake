#set default compiler to use UNIX, compatible with make
#unresolved related issue: when Jan tried opening VS2017-compiled DLL he was missing all 
#kinds of libraries. Does the Unix compiler do better?
set (CMAKE_GENERATOR "Unix Makefiles" CACHE INTERNAL "" FORCE)
