import { useEffect, useRef } from "react";
import igv from "igv";
import { api, ConditionRef } from "./api";

// Embedded igv.js reproducing Paraphase's recommended view:
// group reads by the HP tag, color by the YC tag (native RGB tag).
export function IgvPanel({
  sampleId,
  condition,
  igvGenome,
}: {
  sampleId: string;
  condition: ConditionRef;
  igvGenome: string;
}) {
  const hostRef = useRef<HTMLDivElement>(null);
  const browserRef = useRef<any>(null);

  useEffect(() => {
    let cancelled = false;
    if (!hostRef.current) return;

    const options: any = {
      genome: igvGenome,
      locus: condition.viewer_locus || condition.region,
      tracks: [
        {
          name: `${sampleId} — ${condition.gene} (HP/YC phased)`,
          type: "alignment",
          format: "bam",
          url: api.bamUrl(sampleId),
          indexURL: api.baiUrl(sampleId),
          colorBy: "tag:YC",
          groupBy: "tag:HP",
          showCoverage: true,
          // Squished reads so all haplotypes fit without heavy scrolling.
          displayMode: "SQUISHED",
          // Paraphase realigns the whole locus (tens of kb); render alignments
          // regardless of zoom rather than igv.js's default 30 kb window.
          visibilityWindow: -1,
          height: 460,
        },
        {
          name: `${condition.gene} variants`,
          type: "variant",
          format: "vcf",
          url: api.vcfUrl(sampleId, condition.gene),
          displayMode: "EXPANDED",
          visibilityWindow: -1,
          height: 90,
        },
      ],
    };

    const locus = condition.viewer_locus || condition.region;
    igv.createBrowser(hostRef.current, options).then((b: any) => {
      if (cancelled) {
        igv.removeBrowser(b);
        return;
      }
      browserRef.current = b;
      // Force the intended locus after load — guards against igv auto-fitting
      // to a wider range when the container width isn't settled at creation.
      try {
        b.search(locus);
      } catch {
        /* ignore */
      }
    });

    return () => {
      cancelled = true;
      if (browserRef.current) {
        igv.removeBrowser(browserRef.current);
        browserRef.current = null;
      }
    };
  }, [sampleId, condition.gene, condition.region, igvGenome]);

  return (
    <div>
      <div className="igv-hint">
        Reads grouped by <b>HP</b> (haplotype) and colored by <b>YC</b> tag — the
        Paraphase-recommended view. Blue = uniquely assigned; green/purple = the two
        alleles; gray = unassigned / multi-mapping.
      </div>
      <div className="igv-host" ref={hostRef} />
    </div>
  );
}
