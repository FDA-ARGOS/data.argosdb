import re
import click
import questionary as q
from datetime import datetime
from pathlib import Path
from typing import Set, List, Dict, Tuple, Optional, Union, Iterable, Any
from collections import namedtuple
from questionary.constants import DEFAULT_STYLE

from prompt_toolkit import print_formatted_text, HTML
from prompt_toolkit.validation import Validator, ValidationError
from prompt_toolkit.document import Document

from ..cli import cli
from ..utils.csvv import load_csv, dump_csv

from .new_study import abort_if_ref_name_unused


def load_subjects(payload_dir: Path) -> List[Dict]:
    subjects_csv: Path = payload_dir / 'tables' / 'subjects.csv'
    return load_csv(subjects_csv)


def split_patient_names(text_patient_names: str) -> Iterable[str]:
    pn: str
    for pn in re.split(r'[\r\n,]', text_patient_names):
        pn = pn.strip()
        if not pn:
            continue
        yield pn


def generate_patient_names(
    name_format: str,
    num_patients: int
) -> Iterable[str]:
    if name_format == 'digit':
        for idx in range(num_patients):
            yield 'Patient {}'.format(idx + 1)
    else:  # alphabet
        for idx in range(num_patients):
            yield 'Patient {}'.format(chr(65 + idx))


class ManualPatientNamesValidator(Validator):

    def validate(self, doc: Document) -> None:
        patient_names: Dict[str, str] = {}
        for pn in split_patient_names(doc.text):
            if pn in patient_names:
                raise ValidationError(
                    message=f'Duplicate patient name {pn} supplied'
                )
            pnlower = pn.lower()
            if pnlower in patient_names:
                raise ValidationError(
                    message='Patient name is too similar to {}: {}'
                    .format(patient_names[pnlower], pn)
                )
            patient_names[pnlower] = pn


class NumPatientsValidator(Validator):

    name_format: str

    def __init__(self, name_format: str, *args: Any, **kwargs: Any):
        self.name_format = name_format
        super().__init__(*args, **kwargs)

    def validate(self, doc: Document) -> None:
        text_num: str = doc.text
        if not text_num.isdigit():
            raise ValidationError(
                message=f'A number must be entered, not {text_num}'
            )
        num: int = int(text_num)
        if num < 1:
            raise ValidationError(
                message=f'Number of patients cannot be {num}'
            )
        if self.name_format == 'alphabet' and num > 26:
            raise ValidationError(
                message='# patients cannot be greater 26 when '
                'alphabet format (Patient A, B, ...) was selected'
            )


@cli.command()
@click.argument(
    'payload_dir',
    is_eager=True,
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
@click.option(
    '--ref-name', prompt="Enter refName",
    callback=abort_if_ref_name_unused,
    help='Reference name of article contains data of the subject(s)')
def add_patients(
    payload_dir: str,
    ref_name: str
) -> None:
    tables_dir: Path = Path(payload_dir) / 'tables'

    has_patients: bool = q.confirm(
        f'Does {ref_name} contains (human) patient data? (e.g. age, '
        'infection, virus/plasma sample collection, vaccination, treatment)'
    ).ask()

    if not has_patients:
        return

    name_format: str = q.select(
        'Which format should be used for the patient names?',
        choices=[
            q.Choice('Patient 1, Patient 2, ...', 'digit'),
            q.Choice('Patient A, Patient B, ...', 'alphabet'),
            q.Choice('Enter them one by one', 'manual')
        ]
    ).ask()

    num_patients: int = 0
    patient_names: List[str] = []
    confirm_num_patients: bool = False

    while not confirm_num_patients and name_format == 'manual':
        text_patient_names: str = q.text(
            'Enter the patient names: (separate by commas or line breaks)',
            validate=ManualPatientNamesValidator,
            multiline=True
        ).ask()
        patient_names = list(split_patient_names(text_patient_names))
        num_patients = len(patient_names)
        confirm_num_patients = q.confirm(
            'You have entered n={} patient{}, correct?'.format(
                num_patients,
                's' if num_patients > 1 else ''
            )
        ).ask()

    while not confirm_num_patients and name_format != 'manual':
        text_num_patients: str = q.text(
            f'How many patients does {ref_name} have?',
            validate=NumPatientsValidator(name_format)
        ).ask()
        num_patients = int(text_num_patients)
        patient_names = list(generate_patient_names(name_format, num_patients))
        ellipsis_printed: bool = False
        for idx, pn in enumerate(patient_names):
            if idx < 3 or idx > num_patients - 3:
                print_formatted_text(
                    HTML(f'  <answer>{pn}</answer>'),
                    style=DEFAULT_STYLE
                )
            elif ellipsis_printed is False:
                print_formatted_text(
                    HTML('  <answer>...</answer>'),
                    style=DEFAULT_STYLE
                )
                ellipsis_printed = True
        confirm_num_patients = q.confirm(
            'n={} patient{} generated as above shown. Confirm?'.format(
                num_patients,
                's were' if num_patients > 1 else ' was'
            )
        ).ask()
