import logging
import os
from dataclasses import asdict
from dataclasses import dataclass

from cookiecutter.main import cookiecutter
from sc_toolbox.cli.questionary import sc_toolbox_questionary

log = logging.getLogger(__name__)


@dataclass
class TemplateAttributes:
    template_type: str = ""

    project_name: str = ""
    project_slug: str = ""


class ProjectCreator:
    def __init__(self):
        self.WD = os.path.dirname(__file__)
        self.TEMPLATES_PATH = f"{self.WD}/templates"
        self.template_attributes = TemplateAttributes()

    def create_project(self):
        """
        Prompts and guides the user through a number of possible dataloader choices.
        Prompts the user for required attributes which must be present in the dataloader.
        Finally creates the specific cookiecutter dataloader template.
        """
        self._prompt_project_template()

        switcher = {
            "simple_docker": self._prompt_simple_docker,
        }
        switcher.get(self.template_attributes.template_type)()  # type: ignore

        self._create_dataloader_template()

    def _prompt_project_template(self) -> None:
        """
        Guides the user to select the appropriate dataloader template for his dataset.
        Sets the dataloader_type
        """
        template_type = sc_toolbox_questionary(
            function="select", question="Select a project template:", choices=["Simple Docker"]
        )
        if template_type == "Simple Docker":
            self.template_attributes.template_type = "simple_docker"

    def _prompt_simple_docker(self):
        """
        Prompts the user for all required attributes for a dataloader such as DOI, author, etc.
        """
        self.template_attributes.project_name = sc_toolbox_questionary(
            function="text", question="Project Name", default="Single Cell Analysis"
        )
        self.template_attributes.project_slug = self.template_attributes.project_name.replace(" ", "_")

    def _template_attributes_to_dict(self) -> dict:
        """
        Create a dict from the our Template Structure dataclass
        :return: The dict containing all key-value pairs with non empty values
        """
        return {key: val for key, val in asdict(self.template_attributes).items() if val != ""}

    def _create_dataloader_template(self):
        template_path = f"{self.TEMPLATES_PATH}/{self.template_attributes.template_type}"
        cookiecutter(
            f"{template_path}",
            no_input=True,
            overwrite_if_exists=True,
            extra_context=self._template_attributes_to_dict(),
        )
