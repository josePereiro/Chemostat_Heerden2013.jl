## ------------------------------------------------------------------------
cd(@__DIR__)
# Install unregistered packages
using Pkg
Pkg.activate(Base.current_project(@__DIR__))
try
    pkg"rm Chemostat"
    pkg"rm UtilsJL"
catch; end
pkg"add https://github.com/josePereiro/UtilsJL.git#master"
pkg"add https://github.com/josePereiro/Chemostat#9b2070a"
pkg"instantiate"
pkg"build"
pkg"test Chemostat"