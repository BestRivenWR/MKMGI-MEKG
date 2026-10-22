"""Run the registered five-condition experiment, then report and deployment refit."""
import os,subprocess,sys
from common import CODE,ROOT,config,verify_freeze

if __name__=='__main__':
    verify_freeze(check_sources=True)
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    for condition in config()['conditions80']:
        for fold in range(5):
            subprocess.run([sys.executable,'-B',str(CODE/'train.py'),'--grouping','80',
                            '--fold',str(fold),'--condition',condition],env=env,check=True)
    subprocess.run([sys.executable,'-B',str(CODE/'report.py')],env=env,check=True)
    subprocess.run([sys.executable,'-B',str(CODE/'finalize.py')],env=env,check=True)
