import argparse
from pathlib import Path
from PIL import Image, ImageDraw
import numpy as np
import pandas as pd


def save_grayscale_image(image_path: Path, output_path: Path) -> Image.Image:
    """Save the grayscale version of an image and return it."""
    image = Image.open(image_path).convert("L")
    image.save(output_path)
    print(f"[OK] Grayscale image saved: {output_path}")
    return image


def save_gridded_image(image: Image.Image, grid_rows: int, grid_cols: int, output_path: Path) -> None:
    """Draw a grid overlay on the image and save it."""
    draw = ImageDraw.Draw(image)
    image_width, image_height = image.size
    cell_width = image_width / grid_cols
    cell_height = image_height / grid_rows

    for r in range(1, grid_rows):
        y = int(round(r * cell_height))
        draw.line([(0, y), (image_width, y)], fill="red", width=1)
    for c in range(1, grid_cols):
        x = int(round(c * cell_width))
        draw.line([(x, 0), (x, image_height)], fill="red", width=1)

    image.save(output_path)
    print(f"[OK] Gridded overlay saved: {output_path}")


def extract_pixel_intensities(image: Image.Image, grid_rows: int, grid_cols: int) -> list[list[float]]:
    """
    Extract inverted mean pixel intensities for each grid cell.
    (255 - mean grayscale) scores darker regions higher.
    Returns rows: [row_index(1-based), col_index(1-based), intensity]
    """
    image_width, image_height = image.size
    cell_width = image_width / grid_cols
    cell_height = image_height / grid_rows

    intensities = []
    for r in range(grid_rows):
        for c in range(grid_cols):
            left = int(round(c * cell_width))
            upper = int(round(r * cell_height))
            right = int(round((c + 1) * cell_width))
            lower = int(round((r + 1) * cell_height))

            cell = image.crop((left, upper, right, lower))
            mean_intensity = 255 - float(np.mean(np.array(cell)))
            intensities.append([r + 1, c + 1, mean_intensity])

    return intensities


def parse_img_arg(arg: str) -> tuple[str, Path]:
    """
    Parse --img arguments of the form:
      Label=/path/to/image.tif
    """
    if "=" not in arg:
        raise ValueError(f"Invalid --img '{arg}'. Use Label=/path/to/image.tif")
    label, path = arg.split("=", 1)
    label = label.strip()
    path = Path(path).expanduser().resolve()
    if not label:
        raise ValueError("Image label is empty.")
    if not path.exists():
        raise FileNotFoundError(f"Image path does not exist: {path}")
    return label, path


def safe_colname(label: str) -> str:
    """Make a CSV-safe column name."""
    return (
        label.strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace("(", "")
        .replace(")", "")
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract inverted mean pixel intensities from gridded Typhoon images (any conditions)."
    )
    parser.add_argument(
        "--img",
        action="append",
        required=True,
        help="Condition label + image path in the form Label=/path/to/image.tif . Repeat for multiple conditions.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to save the output CSV file (e.g., intensities.csv)."
    )
    parser.add_argument(
        "--grid", default="16x24", help="Grid size as ROWSxCOLS (default: 16x24)."
    )
    args = parser.parse_args()

    # Parse grid size
    try:
        grid_rows, grid_cols = map(int, args.grid.lower().split("x"))
    except ValueError as e:
        raise SystemExit("Error: --grid must be in the format ROWSxCOLS (e.g., 16x24).") from e

    # Parse output path and ensure output dir exists
    out_csv = Path(args.output).expanduser().resolve()
    out_dir = out_csv.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Parse images
    labels_paths = [parse_img_arg(x) for x in args.img]

    # Process each image -> intensity table keyed by (Row, Column)
    merged = None

    for label, img_path in labels_paths:
        colname = safe_colname(label)

        gray_path = out_dir / f"{colname}_grayscale.png"
        grid_path = out_dir / f"{colname}_gridded.png"

        img_gray = save_grayscale_image(img_path, gray_path)
        save_gridded_image(img_gray.copy(), grid_rows, grid_cols, grid_path)

        intensities = extract_pixel_intensities(img_gray, grid_rows, grid_cols)
        df = pd.DataFrame(intensities, columns=["Row", "Column", colname])

        merged = df if merged is None else pd.merge(merged, df, on=["Row", "Column"], how="inner")

    # Save merged CSV
    merged.to_csv(out_csv, index=False)
    print(f"[OK] Pixel intensities saved: {out_csv}")
    print(f"[INFO] Conditions: {', '.join([safe_colname(l) for l, _ in labels_paths])}")
    print(f"[INFO] Grid: {grid_rows}x{grid_cols}")


if __name__ == "__main__":
    main()
