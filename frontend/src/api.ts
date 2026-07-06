// Thin API client for the FastAPI backend.

export interface SampleSummary {
  sample_id: string;
  build: string;
  genes: string[];
  sample_sex: string | null;
  genome_depth: number | null;
  json_md5: string;
  bam_md5: string;
  paraphase_version: string | null;
  has_bam: boolean;
  warnings: string[];
  indexed_at: string | null;
}

export interface ConditionRef {
  gene: string;
  disease: string;
  summary: string;
  region: string;
  viewer_locus?: string;
}

export interface Evidence {
  label: string;
  value: unknown;
  note?: string | null;
}

export interface Provenance {
  app_version: string;
  paraphase_version: string | null;
  genome_build: string;
  sample_id: string;
  json_md5: string;
  bam_md5: string;
  generated_at: string;
}

export interface Result {
  gene: string;
  disease: string;
  module_version: string;
  headline: string;
  status: string;
  level: string;
  interpretation: string;
  raw: Record<string, unknown>;
  evidence: Evidence[];
  references: string[];
  qc: {
    region_depth_median: number | null;
    region_depth_p80: number | null;
    genome_depth: number | null;
    min_region_depth: number;
    flags: string[];
    pass: boolean;
  };
  caveats: string[];
  experimental?: boolean;
  secondary?: Result | null;
  secondary_tab_label?: string | null;
  provenance?: Provenance;
  condition?: { gene: string; disease: string; region: string; build: string };
}

async function j<T>(url: string, opts?: RequestInit): Promise<T> {
  const r = await fetch(url, opts);
  if (!r.ok) throw new Error(`${r.status} ${await r.text()}`);
  return r.json() as Promise<T>;
}

export const api = {
  health: () => j<{ app_version: string; results_root: string; igv_genome: string }>("/api/health"),
  rescan: () => j<{ indexed: number; samples: string[] }>("/api/rescan", { method: "POST" }),
  samples: () => j<{ samples: SampleSummary[] }>("/api/samples"),
  sampleDetail: (id: string) =>
    j<SampleSummary & { conditions: ConditionRef[] }>(`/api/samples/${id}`),
  result: (id: string, gene: string) => j<Result>(`/api/samples/${id}/${gene}/result`),
  sessionUrl: (id: string, gene: string) => `/api/samples/${id}/${gene}/session.xml`,
  bundleUrl: (id: string, gene: string) => `/api/samples/${id}/${gene}/bundle.zip`,
  rawUrl: (id: string) => `/api/samples/${id}/raw`,
  bamUrl: (id: string) => `/api/files/${id}/bam`,
  baiUrl: (id: string) => `/api/files/${id}/bai`,
  vcfUrl: (id: string, gene: string) => `/api/files/${id}/vcf/${gene}`,
};
