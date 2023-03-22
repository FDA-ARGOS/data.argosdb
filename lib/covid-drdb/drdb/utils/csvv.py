import csv
from pathlib import Path
from itertools import chain
from more_itertools import unique_everseen
from typing import List, Dict, Optional, Union, Any, Iterable

CSVReaderRow = Dict[str, Optional[str]]
CSVWriterRow = Dict[str, Any]


def load_csv(
    file_path: Union[str, Path],
    null_str: str = r'NULL'
) -> List[CSVReaderRow]:
    with open(file_path, encoding='utf-8-sig') as fd:
        row: CSVReaderRow
        rows: List[CSVReaderRow] = []
        for _row in csv.DictReader(fd):
            row = {**_row}
            for key, val in row.items():
                if val == null_str:
                    row[key] = None
            rows.append(row)
        return rows


def load_multiple_csvs(
    csv_dir: Union[str, Path],
    null_str: str = r'NULL'
) -> List[CSVReaderRow]:
    child: Path
    rows: List[CSVReaderRow] = []
    for child in sorted(Path(csv_dir).iterdir()):
        if child.suffix.lower() != '.csv':
            continue
        rows.extend(load_csv(child, null_str))
    return rows


def dump_csv(
    file_path: Union[str, Path],
    records: Iterable[CSVWriterRow],
    headers: List[str] = [],
    BOM: bool = False,
    null_str: str = r'NULL'
) -> None:
    encoding: str
    writer: csv.DictWriter
    key: str
    row: CSVWriterRow
    _records: List[CSVWriterRow] = list(records)
    if not _records:
        return
    if not headers:
        headers = list(unique_everseen(
            chain(*[r.keys() for r in _records])
        ))

    if BOM:
        encoding = 'utf-8-sig'
    else:
        encoding = 'utf-8'

    with open(file_path, 'w', encoding=encoding) as fd:
        writer = csv.DictWriter(
            fd,
            fieldnames=headers,
            restval=null_str,
            extrasaction='ignore')
        writer.writeheader()
        for row in _records:
            for key, val in row.items():
                if val is None:
                    row[key] = null_str
            writer.writerow(row)
