import re
import logging
from typing import Dict, Any

from scm.plams.core.logging import CSVFormatter
from scm.plams.core.basejob import Job
from scm.plams.core.enums import JobStatus
from scm.plams.core.errors import ResultsError

__all__ = ["JobCSVFormatter"]


class JobCSVFormatter(CSVFormatter):
    """
    Formatter which creates comma-separated log lines from a ``Job`` log record.
    """

    def format(self, record: logging.LogRecord) -> str:
        if isinstance(record.msg, Job):
            record.msg = self._format_job(record.msg)
        return super().format(record)

    @staticmethod
    def _format_job(job: Job) -> Dict[str, Any]:
        message = {
            "job_base_name": re.sub(r"\.\d+$", "", job.name),
            "job_name": job.name,
            "job_status": job.status,
            "job_path": "",
            "job_ok": "",
            "job_check": "",
            "job_get_errormsg": "",
            "job_timeline": str.join(" -> ", [f"{dt.strftime('%Y-%m-%d %H:%M:%S')} {s}" for dt, s in job.status_log]),
            "job_parent_name": "",
            "job_parent_path": "",
        }

        if job.status not in [JobStatus.CREATED, JobStatus.STARTED]:
            message.update({"job_path": job.path})

            if job.status not in [JobStatus.REGISTERED, JobStatus.RUNNING, JobStatus.DELETED]:
                # Avoid race condition when accessing results for deleted job
                try:
                    message.update({"job_ok": job.ok()})
                    message.update({"job_check": job.check()})
                    message.update({"job_get_errormsg": job.get_errormsg()})
                except ResultsError:
                    pass

        if job.parent:
            message["job_parent_name"] = job.parent.name
            message["job_parent_path"] = job.parent.path

        return message
