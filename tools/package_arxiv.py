"""Package the manuscript sources and a separate reading copy for arXiv."""

from pathlib import Path
from shutil import copy2
from zipfile import ZIP_DEFLATED, ZipFile


def main():
    root = Path(__file__).resolve().parents[1]
    paper = root / "paper"
    pdf = paper / "main.pdf"
    if not pdf.is_file():
        raise SystemExit("Compile the manuscript first with make paper.")

    output = root / "output" / "pdf"
    output.mkdir(parents=True, exist_ok=True)
    archive = output / "order15_arxiv_source.zip"
    sources = sorted(paper.rglob("*.tex"))
    with ZipFile(archive, "w", compression=ZIP_DEFLATED) as package:
        for source in sources:
            package.write(source, source.relative_to(paper).as_posix())
    copy2(pdf, output / "order15_arxiv.pdf")
    print(f"Created {archive.relative_to(root)} ({len(sources)} source files)")
    print("Created output/pdf/order15_arxiv.pdf (separate reading copy)")


if __name__ == "__main__":
    main()
