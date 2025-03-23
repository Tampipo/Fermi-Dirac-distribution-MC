from Args import SetupArgs, OrganizeArgs

if __name__ == '__main__':
    args = SetupArgs()
    job_args = OrganizeArgs(args)
    print(job_args)
    if job_args['job'] == 'part_scan':
        from part_scan import main
        main(job_args)
    elif job_args['job'] == 'temp_scan':
        from temp_scan import main
        main(job_args)
    else:    
        from phase_space_scan import main
        main(job_args)