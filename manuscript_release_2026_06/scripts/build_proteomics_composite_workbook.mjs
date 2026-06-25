import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(scriptDir, "..");
const tag = process.argv.includes("--tag")
  ? process.argv[process.argv.indexOf("--tag") + 1]
  : "current";

const compositeDir = path.join(root, "evidence_tables", "composite_proteomics_priority");
const litDir = path.join(root, "evidence_tables", "literature_proteomics_family_mapping");
const outputDir = path.join(root, "outputs", `composite_proteomics_priority_${tag}`);
await fs.mkdir(outputDir, { recursive: true });

const workbook = Workbook.create();

function columnName(index) {
  let n = index + 1;
  let out = "";
  while (n > 0) {
    const rem = (n - 1) % 26;
    out = String.fromCharCode(65 + rem) + out;
    n = Math.floor((n - 1) / 26);
  }
  return out;
}

function coerceCell(value) {
  if (value === "") return null;
  if (value === "True") return true;
  if (value === "False") return false;
  if (/^-?\d+(?:\.\d+)?(?:e[+-]?\d+)?$/i.test(value)) return Number(value);
  return value;
}

function parseCsv(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let inQuotes = false;
  for (let i = 0; i < text.length; i += 1) {
    const ch = text[i];
    const next = text[i + 1];
    if (inQuotes) {
      if (ch === '"' && next === '"') {
        cell += '"';
        i += 1;
      } else if (ch === '"') {
        inQuotes = false;
      } else {
        cell += ch;
      }
    } else if (ch === '"') {
      inQuotes = true;
    } else if (ch === ",") {
      row.push(cell);
      cell = "";
    } else if (ch === "\n") {
      row.push(cell);
      rows.push(row);
      row = [];
      cell = "";
    } else if (ch !== "\r") {
      cell += ch;
    }
  }
  if (cell.length > 0 || row.length > 0) {
    row.push(cell);
    rows.push(row);
  }
  const width = Math.max(0, ...rows.map((r) => r.length));
  return rows.map((r) => {
    const out = r.map(coerceCell);
    while (out.length < width) out.push(null);
    return out;
  });
}

function writeRows(sheet, rows) {
  if (!rows.length || !rows[0].length) return;
  const last = `${columnName(rows[0].length - 1)}${rows.length}`;
  sheet.getRange(`A1:${last}`).values = rows;
}

function styleSheet(sheet) {
  sheet.showGridLines = false;
  try {
    sheet.freezePanes.freezeRows(1);
  } catch {}
  const used = sheet.getUsedRange(true);
  if (!used) return;
  used.format.font = { name: "Aptos", size: 10 };
  used.format.wrapText = false;
  try {
    used.format.autofitColumns();
    used.format.autofitRows();
  } catch {}
  try {
    used.getRow(0).format = {
      fill: "#1F4E79",
      font: { bold: true, color: "#FFFFFF" },
      wrapText: true,
      verticalAlignment: "top",
    };
    used.format.borders = {
      insideHorizontal: { style: "thin", color: "#E5E7EB" },
      insideVertical: { style: "thin", color: "#F3F4F6" },
      bottom: { style: "thin", color: "#CBD5E1" },
    };
  } catch {}
}

async function importCsvSheet(filePath, sheetName) {
  const text = await fs.readFile(filePath, "utf8");
  const rows = parseCsv(text);
  const sheet = workbook.worksheets.add(sheetName);
  writeRows(sheet, rows);
  styleSheet(sheet);
  return sheet;
}

const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
readme.getRange("A1:H1").merge();
readme.getRange("A1").values = [["Composite proteomics-prioritized family evidence"]];
readme.getRange("A1").format = {
  fill: "#0F766E",
  font: { bold: true, color: "#FFFFFF", size: 14 },
};
readme.getRange("A3:B13").values = [
  ["Generated tag", tag],
  ["Purpose", "Join current Model-Supported / Highly Pure protein-family inventories to external literature proteomics filtered for active-phase-up nitrogen-fixation evidence, related protein-family atlas concordance, condensate ranking, and morphotype breadth."],
  ["Primary score", "scope_adjusted_story_score = accessory_story_score plus a Highly Pure related-atlas scope adjustment where collaborator-atlas evidence was not applicable."],
  ["Interpretation", "The score is a prioritization/down-filter for biologically story-rich families, not a formal validation statistic."],
  ["External proteomics", "Scores independent accession-mapped source rows only when they increase in the organism-appropriate N-fix-active phase/cell. Broad response, active-phase-down, transition, and non-target phase rows remain visible as metadata but do not add points."],
  ["Excluded layers", "Local Cyanothece proteomics and prior FOX probability/status layers are excluded from both scoring and visible composite output."],
  ["Related atlas", "Captures convergence with the collaborator/independent protein-family atlas at product/function level."],
  ["Best first tab", "TopFamilies"],
  ["Full joined table", "AllComposite"],
  ["High-purity markers", "HighPurityMarkers lists evidence-supported Highly Pure families that pass the N-fix-response/score filter."],
  ["Raw-ish evidence tabs", "LiteratureNfixResponse, ExternalProteinRows"],
];
readme.getRange("A3:A13").format = { font: { bold: true }, fill: "#E0F2FE" };
readme.getRange("B3:B13").format = { wrapText: true };
readme.getRange("A:A").format.columnWidth = 24;
readme.getRange("B:B").format.columnWidth = 110;

await importCsvSheet(path.join(compositeDir, "summary_metrics.csv"), "Summary");
await importCsvSheet(path.join(compositeDir, "top_story_families.csv"), "TopFamilies");
await importCsvSheet(path.join(compositeDir, "high_purity_marker_families.csv"), "HighPurityMarkers");
await importCsvSheet(path.join(compositeDir, "module_story_summary.csv"), "ModuleSummary");
await importCsvSheet(path.join(compositeDir, "component_definitions.csv"), "ScoreRules");
await importCsvSheet(path.join(compositeDir, "family_composite_all.csv"), "AllComposite");
await importCsvSheet(path.join(compositeDir, "external_literature_nfix_response_evidence.csv"), "LiteratureNfixResponse");
await importCsvSheet(path.join(compositeDir, "external_literature_protein_rows_nfix_classified.csv"), "ExternalProteinRows");
await importCsvSheet(path.join(litDir, "mapping_summary.csv"), "MappingSummary");

for (const sheetName of ["TopFamilies", "HighPurityMarkers", "AllComposite", "LiteratureNfixResponse", "ExternalProteinRows"]) {
  try {
    workbook.worksheets.getItem(sheetName).freezePanes.freezeColumns(2);
  } catch {}
}

for (const sheetName of ["Summary", "TopFamilies", "HighPurityMarkers", "ModuleSummary", "ScoreRules"]) {
  const preview = await workbook.render({
    sheetName,
    range: "A1:J25",
    scale: 1,
    format: "png",
  });
  await fs.writeFile(
    path.join(outputDir, `${sheetName}.preview.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const summaryInspect = await workbook.inspect({
  kind: "table",
  range: "Summary!A1:B20",
  include: "values,formulas",
  tableMaxRows: 20,
  tableMaxCols: 3,
  maxChars: 4000,
});
console.log(summaryInspect.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

const output = await SpreadsheetFile.exportXlsx(workbook);
const outputPath = path.join(outputDir, "proteomics_composite_family_evidence.xlsx");
await output.save(outputPath);
console.log(`SAVED ${outputPath}`);
