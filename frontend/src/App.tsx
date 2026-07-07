import { useEffect, useMemo, useState } from "react";
import {
  api,
  ConditionRef,
  Result,
  SampleSummary,
} from "./api";
import { IgvPanel } from "./IgvPanel";

type Detail = SampleSummary & { conditions: ConditionRef[] };
type Tab = "results" | "carrier" | "viewer" | "download" | "raw";

// Global condition filter. "all" shows every condition a sample has; picking a
// specific condition restricts the sample list and the selected sample to just
// that condition.
const CONDITION_OPTIONS: { key: string; label: string }[] = [
  { key: "all", label: "All" },
  { key: "smn1", label: "SMA" },
  { key: "f8", label: "Hemophilia A" },
];

export function App() {
  const [samples, setSamples] = useState<SampleSummary[]>([]);
  const [filter, setFilter] = useState("");
  const [conditionFilter, setConditionFilter] = useState("all");
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [detail, setDetail] = useState<Detail | null>(null);
  const [gene, setGene] = useState<string | null>(null);
  const [result, setResult] = useState<Result | null>(null);
  const [tab, setTab] = useState<Tab>("results");
  const [igvGenome, setIgvGenome] = useState("hg38");
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    api.health().then((h) => setIgvGenome(h.igv_genome)).catch(() => {});
    loadSamples();
  }, []);

  function loadSamples() {
    api.samples().then((r) => setSamples(r.samples)).catch((e) => setError(String(e)));
  }

  async function rescan() {
    setError(null);
    try {
      await api.rescan();
      loadSamples();
    } catch (e) {
      setError(String(e));
    }
  }

  useEffect(() => {
    if (!selectedId) return;
    setDetail(null);
    setResult(null);
    setGene(null);
    api.sampleDetail(selectedId).then(setDetail);
  }, [selectedId]);

  // Conditions visible for the current sample, honoring the global filter.
  const visibleConditions = useMemo(() => {
    if (!detail) return [];
    return detail.conditions.filter(
      (c) => conditionFilter === "all" || c.gene === conditionFilter
    );
  }, [detail, conditionFilter]);

  // Keep the selected gene valid as the sample / filter changes.
  useEffect(() => {
    if (!detail) return;
    if (visibleConditions.length === 0) {
      setGene(null);
    } else if (!visibleConditions.some((c) => c.gene === gene)) {
      setGene(visibleConditions[0].gene);
    }
  }, [detail, visibleConditions]);

  useEffect(() => {
    if (!selectedId || !gene) return;
    setResult(null);
    setTab("results");
    api.result(selectedId, gene).then(setResult).catch((e) => setError(String(e)));
  }, [selectedId, gene]);

  const filtered = useMemo(
    () =>
      samples.filter(
        (s) =>
          s.sample_id.toLowerCase().includes(filter.toLowerCase()) &&
          (conditionFilter === "all" || s.genes.includes(conditionFilter))
      ),
    [samples, filter, conditionFilter]
  );

  // If the active sample no longer matches the condition filter, clear it.
  useEffect(() => {
    if (selectedId && !filtered.some((s) => s.sample_id === selectedId)) {
      setSelectedId(null);
      setDetail(null);
    }
  }, [filtered, selectedId]);

  const condition = detail?.conditions.find((c) => c.gene === gene) || null;

  return (
    <div className="app">
      <aside className="sidebar">
        <div className="sidebar-header">
          <h1>Paraphase Clinical Review</h1>
          <div className="sub">Long-read segmental-duplication caller · genome {igvGenome}</div>
        </div>
        <div className="toolbar">
          <input
            placeholder="Filter samples…"
            value={filter}
            onChange={(e) => setFilter(e.target.value)}
          />
          <div className="segmented" style={{ marginTop: 10 }}>
            {CONDITION_OPTIONS.map((o) => (
              <button
                key={o.key}
                className={"seg-btn" + (o.key === conditionFilter ? " active" : "")}
                onClick={() => setConditionFilter(o.key)}
              >
                {o.label}
              </button>
            ))}
          </div>
          <div style={{ marginTop: 10 }}>
            <button className="btn secondary" onClick={rescan}>Rescan folder</button>
          </div>
        </div>
        <div className="sample-list">
          {filtered.map((s) => (
            <div
              key={s.sample_id}
              className={"sample-item" + (s.sample_id === selectedId ? " active" : "")}
              onClick={() => setSelectedId(s.sample_id)}
            >
              <div className="sid">{s.sample_id}</div>
              <div className="meta">
                {s.genes.map((g) => (
                  <span className="chip" key={g}>{g}</span>
                ))}
                build {s.build}{s.sample_sex ? ` · ${s.sample_sex}` : ""}
              </div>
            </div>
          ))}
          {filtered.length === 0 && <div className="placeholder">No samples indexed.</div>}
        </div>
      </aside>

      <main className="main">
        {error && <div className="banner warn">{error}</div>}
        {!selectedId && <div className="placeholder">Select a sample to begin.</div>}

        {detail && (
          <>
            <div className="conditions">
              {visibleConditions.map((c) => (
                <button
                  key={c.gene}
                  className={"cond-btn" + (c.gene === gene ? " active" : "")}
                  onClick={() => setGene(c.gene)}
                >
                  {c.disease}
                </button>
              ))}
              {visibleConditions.length === 0 && (
                <div className="placeholder">No interpreted conditions for this sample.</div>
              )}
            </div>

            {condition && (
              <div className="tabs">
                {([
                  { key: "results", label: result?.secondary ? "Status" : "Results" },
                  ...(result?.secondary
                    ? [{ key: "carrier", label: result.secondary_tab_label || "Silent carrier (experimental)" }]
                    : []),
                  { key: "viewer", label: "Visualization" },
                  { key: "download", label: "Download" },
                  { key: "raw", label: "Raw JSON" },
                ] as { key: Tab; label: string }[]).map((t) => (
                  <div
                    key={t.key}
                    className={"tab" + (t.key === tab ? " active" : "")}
                    onClick={() => setTab(t.key)}
                  >
                    {t.label}
                  </div>
                ))}
              </div>
            )}

            {tab === "results" && result && <ResultView result={result} />}
            {tab === "carrier" && result?.secondary && (
              <ResultView result={result.secondary} />
            )}
            {tab === "viewer" && condition && detail.has_bam && (
              <IgvPanel sampleId={detail.sample_id} condition={condition} igvGenome={igvGenome} />
            )}
            {tab === "viewer" && !detail.has_bam && (
              <div className="banner warn">No BAM available for this sample.</div>
            )}
            {tab === "download" && condition && (
              <DownloadTab sampleId={detail.sample_id} condition={condition} />
            )}
            {tab === "raw" && (
              <pre className="raw">{/* loaded lazily below */}<RawJson id={detail.sample_id} /></pre>
            )}
          </>
        )}
      </main>
    </div>
  );
}

function ResultView({ result }: { result: Result }) {
  return (
    <div>
      {result.experimental && (
        <div className="banner warn">
          <b>Experimental — research use only.</b> This silent-carrier (2+0)
          assessment is not validated for clinical reporting.
        </div>
      )}
      <div className={"headline-card " + result.level}>
        <h2>{result.headline}</h2>
        <p>{result.interpretation}</p>
      </div>

      <div className="card">
        <h3>Quality control</h3>
        <p>
          Status:{" "}
          {result.qc.pass ? <span className="qc-pass">PASS</span> : <span className="qc-fail">REVIEW</span>}
          {"  ·  "}region depth (median): {String(result.qc.region_depth_median)}x
          {"  ·  "}genome depth: {String(result.qc.genome_depth)}x
        </p>
        {result.qc.flags.map((f, i) => (
          <div className="banner warn" key={i}>{f}</div>
        ))}
      </div>

      <div className="card">
        <h3>Evidence &amp; raw calls</h3>
        <table className="kv">
          <tbody>
            {result.evidence.map((e, i) => (
              <tr key={i}>
                <td className="k">{e.label}</td>
                <td>
                  <span className="mono">
                    {typeof e.value === "object" ? JSON.stringify(e.value) : String(e.value)}
                  </span>
                  {e.note && <div className="note">{e.note}</div>}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {result.caveats.length > 0 && (
        <div className="card">
          <h3>Caveats</h3>
          {result.caveats.map((c, i) => (
            <p className="caveat" key={i}>• {c}</p>
          ))}
        </div>
      )}

      <div className="card">
        <h3>References</h3>
        <ul className="refs">
          {result.references.map((r, i) => (
            <li key={i}>{r}</li>
          ))}
        </ul>
        {result.provenance && (
          <div className="provbar">
            Interpreter <b>{result.gene}</b> v{result.module_version} ·
            Paraphase {result.provenance.paraphase_version || "?"} ·
            app v{result.provenance.app_version} ·
            build {result.provenance.genome_build} ·
            JSON md5 {result.provenance.json_md5.slice(0, 12)} ·
            BAM md5 {result.provenance.bam_md5?.slice(0, 12)} ·
            generated {result.provenance.generated_at}
          </div>
        )}
      </div>
    </div>
  );
}

function DownloadTab({ sampleId, condition }: { sampleId: string; condition: ConditionRef }) {
  return (
    <div className="card">
      <h3>Download — IGV assets</h3>
      <p className="caveat">
        The IGV desktop session reproduces the Paraphase-recommended view (group by
        HP, color by YC) at {condition.region}. The bundle additionally packs the
        BAM, index and VCF with a portable (relative-path) session so it opens on any
        machine.
      </p>
      <div style={{ marginTop: 12 }}>
        <a className="btn" href={api.sessionUrl(sampleId, condition.gene)}>
          IGV session (.xml)
        </a>
        <a className="btn secondary" href={api.bundleUrl(sampleId, condition.gene)}>
          Portable bundle (.zip)
        </a>
        <a className="btn secondary" href={api.rawUrl(sampleId)} target="_blank">
          Raw Paraphase JSON
        </a>
      </div>
    </div>
  );
}

function RawJson({ id }: { id: string }) {
  const [text, setText] = useState("Loading…");
  useEffect(() => {
    fetch(api.rawUrl(id))
      .then((r) => r.json())
      .then((d) => setText(JSON.stringify(d, null, 2)))
      .catch((e) => setText(String(e)));
  }, [id]);
  return <>{text}</>;
}
