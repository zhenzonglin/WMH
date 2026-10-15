"""Synthetic publication-table contracts; no clinical data are embedded."""
from __future__ import annotations

import csv
import warnings
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("tableone")
pytest.importorskip("docx")
pytest.importorskip("lxml")

from lxml import etree

from wmh_hcy.recovery_publication.tables import (
    DEMO_FOOTNOTE,
    PRIMARY,
    _baseline_content,
    build_baseline_tables,
    export_table,
)


def _cohort() -> pd.DataFrame:
    return pd.DataFrame({
        "patient_id": [f"PRIVATE-PATIENT-{i}" for i in range(8)],
        "state60": [0, 0, 0, 1, 1, 2, 2, np.nan],
        "age": [50, 60, 70, 55, 65, 66, 76, 1000],
        "wmh_ml": [10, 12, 14, 15, 17, 20, 22, 1000],
        "gm119_ml": [550, 560, 570, 530, 540, 510, 520, 1000],
        "icv_ml": [1300] * 8, "lesion_ml": [4, 6, 8, 3, 5, 6, 8, 1000],
        "nihss": [1, 2, 3, 1, 2, 3, 4, 1000],
        "sex": [1, 2, 1, 2, 1, 2, 1, 1],
        "smoking": [1, 2, 3, 4, 1, 2, 3, 4],
        "drinking": [4, 3, 2, 1, 4, 3, 2, 1],
        "hypertension": [1, 2, np.nan, 2, 2, 2, 1, 1],
        "diabetes": [1, 1, 2, 2, 1, 1, 2, 2],
        "prior_stroke": [1, 1, 1, 2, 2, 1, 1, 1],
        "education": [1, 2, np.nan, 3, 4, 5, 2, 1],
        "pre_mrs": [0, 1, 2, 3, 4, 5, 0, 0],
        "toast": [1, 2, 3, 4, 5, 1, 2, 1],
        "mrs3": [0, 1, 2, 0, 1, 2, 0, 0],
    })


def test_tableone_baseline_uses_known_states_nonmissing_denominators():
    columns, rows, notes = _baseline_content(_cohort(), list(PRIMARY))
    assert columns == ["Characteristic", "All known\n(n=7)", "Independent\n(n=3)",
                       "Dependent\n(n=2)", "Dead\n(n=2)"]
    hypertension_index = next(i for i, row in enumerate(rows) if row[0] == "Hypertension")
    assert rows[hypertension_index + 1] == ["  No", "2 (33.3)", "1 (50.0)", "0 (0.0)", "1 (50.0)"]
    assert rows[hypertension_index + 2] == ["  Yes", "4 (66.7)", "1 (50.0)", "2 (100.0)", "1 (50.0)"]
    assert rows[hypertension_index + 3] == ["  Missing, n", "1", "1", "0", "0"]
    age_row = next(row for row in rows if row[0] == "Age, years")
    assert age_row[2] == "60.0 [55.0,65.0]"
    assert not any("1000" in cell or "PATIENT" in cell for row in rows for cell in row)
    assert "nonmissing observations" in " ".join(notes)
    assert "No baseline significance tests" in " ".join(notes)


def test_categories_are_named_without_guessing_education():
    cohort = _cohort()
    _, rows, notes = _baseline_content(cohort, ["sex", "education", "toast"])
    descriptions = {row[0] for row in rows}
    assert "  Male" in descriptions and "  Female" in descriptions
    assert "  Large artery atherosclerosis" in descriptions
    assert "  Education category 1" in descriptions
    assert any("pending verification" in note for note in notes)
    labels = {str(i): f"Verified label {i}" for i in range(1, 6)}
    _, rows, notes = _baseline_content(cohort, ["education"], labels)
    assert "  Verified label 1" in {row[0] for row in rows}
    assert not any("pending verification" in note for note in notes)
    with pytest.raises(ValueError, match="Unmapped category"):
        _baseline_content(cohort, ["education"], {"1": "Incomplete codebook"})


@pytest.mark.parametrize("column,value,message", [
    ("sex", 9, "Unmapped category"), ("age", -1, "Negative values"),
    ("state60", 6, "state60 must use"), ("wmh_ml", np.inf, "infinite values"),
])
def test_invalid_codes_and_measurements_are_not_silently_relabelled(column, value, message):
    cohort = _cohort()
    cohort[column] = cohort[column].astype(float)
    cohort.loc[0, column] = value
    with pytest.raises(ValueError, match=message):
        _baseline_content(cohort, list(PRIMARY))


def test_missing_required_column_blocks_baseline():
    with pytest.raises(ValueError, match="required columns: age"):
        _baseline_content(_cohort().drop(columns="age"), list(PRIMARY))


def test_word_table_has_only_three_rules_and_repeating_header(tmp_path: Path):
    outputs = export_table(tmp_path, "example", "Table 1 Example", ["Characteristic", "Overall"],
                           [["Age, years", "61 [54, 70]"], ["Female", "42 (42.0)"]], ["Values are aggregates."])
    assert {path.suffix for path in outputs} == {".docx", ".csv", ".tex", ".html"}
    with ZipFile(tmp_path / "example.docx") as archive:
        xml = etree.fromstring(archive.read("word/document.xml"))
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    borders = xml.xpath("//w:tblPr/w:tblBorders/*", namespaces=ns)
    actual = {etree.QName(edge).localname: edge.get(f"{{{ns['w']}}}val") for edge in borders}
    assert actual == {"top": "single", "left": "nil", "bottom": "single", "right": "nil",
                      "insideH": "nil", "insideV": "nil"}
    assert len(xml.xpath("//w:tr[1]/w:trPr/w:tblHeader", namespaces=ns)) == 1
    assert len(xml.xpath("//w:tr/w:trPr/w:cantSplit", namespaces=ns)) == 3
    cell_rules = xml.xpath("//w:tcBorders/*", namespaces=ns)
    assert len(cell_rules) == 2
    assert all(etree.QName(rule).localname == "bottom" for rule in cell_rules)
    assert not xml.xpath("//w:shd", namespaces=ns)
    tex = (tmp_path / "example.tex").read_text(encoding="utf-8")
    assert "booktabs" in tex and r"\toprule" in tex and r"\midrule" in tex and r"\bottomrule" in tex
    assert "|" not in tex
    with (tmp_path / "example.csv").open(encoding="utf-8-sig", newline="") as handle:
        assert list(csv.reader(handle))[1] == ["Age, years", "61 [54, 70]"]


def test_baseline_outputs_contain_no_patient_identifiers(tmp_path: Path):
    cohort = _cohort()
    before = cohort.copy(deep=True)
    outputs = build_baseline_tables(cohort, tmp_path)
    pd.testing.assert_frame_equal(cohort, before)
    assert len(outputs) == 8
    for output in outputs:
        if output.suffix == ".docx":
            with ZipFile(output) as archive:
                contents = archive.read("word/document.xml").decode("utf-8")
        else:
            contents = output.read_text(encoding="utf-8-sig")
        assert "patient_id" not in contents
        assert "PRIVATE-PATIENT" not in contents


def test_table_stem_and_ragged_rows_are_rejected(tmp_path: Path):
    with pytest.raises(ValueError, match="Table stem"):
        export_table(tmp_path, "../escape", "Title", ["A"], [["1"]], [])
    with pytest.raises(ValueError, match="same number"):
        export_table(tmp_path, "valid", "Title", ["A", "B"], [["1"]], [])


def test_supplement_invalid_code_blocks_all_baseline_exports(tmp_path: Path):
    cohort = _cohort()
    cohort.loc[0, "education"] = 9
    with pytest.raises(ValueError, match="Unmapped category codes in baseline column education"):
        build_baseline_tables(cohort, tmp_path)
    assert not list(tmp_path.iterdir())


def test_absent_categories_keep_valid_order_without_invented_rows():
    cohort = _cohort().loc[lambda frame: frame.pre_mrs.ne(5)]
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter("always")
        _, rows, _ = _baseline_content(cohort, ["pre_mrs"])
    assert not any("Order value not found" in str(warning.message) for warning in observed)
    assert [row[0] for row in rows] == ["Pre-stroke mRS score", "  0", "  1", "  2", "  3", "  4"]


def test_demo_watermark_is_present_in_every_table_format(tmp_path: Path):
    outputs = build_baseline_tables(_cohort(), tmp_path, demo=True)
    for output in outputs:
        if output.suffix == ".docx":
            with ZipFile(output) as archive:
                contents = archive.read("word/document.xml").decode("utf-8")
        else:
            contents = output.read_text(encoding="utf-8-sig")
        assert "Synthetic style preview" in contents
        assert DEMO_FOOTNOTE in contents
    with (tmp_path / "table1.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    assert all(len(row) == 5 for row in rows)
    assert rows[-1][0].startswith("# Footnote:")
