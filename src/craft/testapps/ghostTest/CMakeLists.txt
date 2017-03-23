add_executable(ghostTest.bin ghostTest.cpp )
add_executable(lanczos.bin lanczos.cpp imtql1.f imtql1f.f pythag.f pythagf.f)
set_target_properties(lanczos.bin PROPERTIES LINKER_LANGUAGE C)


