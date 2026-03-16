import subprocess


class SubmitJobsSlurm(object):
    # init
    def __init__(self, this_logger):
        self.logger = this_logger
        self.output = ""
        self.error = ""
        self.chdir = ""
        self.name = ""
        self.log_jobid = ""
        self.jobid = ""
        self.dependencies = ""

    def set_details(self, name, output, error, chdir):
        self.name = name
        self.output = output
        self.error = error
        self.chdir = chdir

    def set_jobid(self, this_job_id):
        self.jobid = this_job_id

    def set_log_jobid(self, this_job_id_logfile):
        self.log_jobid = this_job_id_logfile

    def set_dependencies(self, this_dependencies):
        self.dependencies = f'--dependency={this_dependencies}'

    def make_submit_job(self, script):
        # with job id in stdout
        job_starter = (f'sbatch {self.dependencies}  --job-name {self.name}  --output {self.output} '
                       f'  --error {self.error}  --chdir {self.chdir} {script}')
        self.logger.info(f'CMD: {job_starter}')
        return job_starter

    # slurm based
    def check_job(self):
        os_call = subprocess.run(f'squeue --long | grep "{self.jobid}"', shell=True, capture_output=True, text=True)
        my_job_status = [info for info in os_call.stdout.split(" ") if info != ""]
        return ",".join(my_job_status)

    def is_job_running(self):
        return self.check_job != ""
