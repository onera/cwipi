set (CWP_DEF_VERSION "1.3.0-a0.untagged")
#set (CWP_DEF_VERSION "1.0.0-a0.untagged")
string(REPLACE "-" ";" CWP_DEF_VERSION_LIST ${CWP_DEF_VERSION})
list(GET CWP_DEF_VERSION_LIST 0 CWP_DEF_VERSION_MAJOR)
