add_executable(dbsr_prep3 mod_dbsr_prep.f90       dbsr_prep.f90
     add_sub_orbital.f90     assign_index.f90      check_bsr_name.f90
     dbsr_prep_inf.f90       read_arg.f90          sub_check_orb.f90
     sub_phys_orb.f90        sub_phys_pert.f90     write_targ_example.f90)
target_link_libraries(dbsr_prep3 dbs zconfjj zcom)
