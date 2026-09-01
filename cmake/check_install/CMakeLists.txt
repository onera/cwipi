cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.16)

project(fake_install LANGUAGES C)

find_package(cwp
             REQUIRED COMPONENTS cwp cwpf
             HINTS $ENV{CWIPI_ROOT}
             )