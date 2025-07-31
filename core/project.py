import os
from typing import Optional, Union, Dict, Literal
from pathlib import Path
import shutil
import json

from scm.plams.core.errors import ProjectError
from scm.plams.core.functions import log


class Project:

    def __init__(
        self,
        name: str,
        description: Optional[str] = None,
        parent_dir: Optional[Union[str, os.PathLike]] = None,
        mode: Literal["x", "r", "r+", "w"] = "x",
    ):
        if not name:
            raise ValueError("Value for 'name' must be specified, it cannot be an empty string or 'None'")
        self._name = name
        self._description = description
        self._path = Path(parent_dir if parent_dir else os.getcwd()) / name
        self._mode = mode

        if self._mode == "x":
            if self._path.exists():
                raise ProjectError(
                    f"Project '{name}' already exists in the directory '{parent_dir}'. Modify the location to create a new project, or change the access mode to use an existing project."
                )
            self.path.mkdir(parents=True)
            self._save()
            log(f"Created project '{name}' in the directory '{parent_dir}'")
        elif self._mode == "r":
            raise NotImplementedError("readonly mode not implemented")
        elif self._mode == "r+":
            raise NotImplementedError("readwrite mode not implemented")
        elif mode == "w":
            if self.path.exists():
                shutil.rmtree(self.path)
                log(f"Deleted existing project '{name}' in the directory '{parent_dir}'")
            self.path.mkdir(parents=True)
            self._save()
            log(f"Created project '{name}' in the directory '{parent_dir}'")
        else:
            raise ValueError(f"Invalid access mode '{mode}', must be one of: 'x' (create), 'r' (read), 'r+' (readwrite) or 'w' (write)")

    @classmethod
    def create(cls, name: str, description: Optional[str] = None, parent_dir: Optional[Union[str, os.PathLike]] = None) -> "Project":
        """
        Create a new project.
        """
        return cls(name, description=description, parent_dir=parent_dir, mode="x")

    @classmethod
    def load(cls, name: str, parent_dir: Optional[Union[str, os.PathLike]] = None):
        raise NotImplementedError("load not implemented")

    # @classmethod
    # def delete(cls, project: "Project"):


    @property
    def name(self) -> str:
        """
        Name of the project. This is also the name of the project directory.

        :return: name of the project
        """
        return self._name

    @property
    def description(self) -> Optional[str]:
        """
        Description of the project and its contents.

        :return: description of the project
        """
        return self._description

    @description.setter
    def description(self, value: Optional[str]):
        self._description = value
        self._save()

    @property
    def path(self) -> Path:
        """
        Absolute path for the project directory.

        :return: full path of the project
        """
        return self._path.resolve()

    @property
    def _metadata_file(self) -> Path:
        return self.path / ".project.json"

    @property
    def _metadata(self) -> Dict:
        return {
            "name": self.name,
            "description": self.description,
        }

    def _save(self):
        metadata = self._metadata
        with open(self._metadata_file, "w") as f:
            json.dump(metadata, f, indent=4)
