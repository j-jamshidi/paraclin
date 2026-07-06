import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// The dev server proxies /api to the FastAPI backend so the browser can
// range-request BAM/VCF from the same origin (required by igv.js). Because the
// proxy runs on the dev-server host, LAN clients that reach this dev server can
// use the backend transparently — no client-side URL changes needed.
//
// Env vars:
//   PARACLIN_HOST     interface to bind (default "localhost"; `npm run dev:lan`
//                     sets 0.0.0.0 via the --host flag so the app is reachable on
//                     your local network).
//   PARACLIN_PORT     dev-server port (default 5199).
//   PARACLIN_BACKEND  backend URL the /api proxy targets (default
//                     http://127.0.0.1:8077). Set this if the backend runs on a
//                     different host/port.
export default defineConfig({
  plugins: [react()],
  server: {
    host: process.env.PARACLIN_HOST || "localhost",
    port: Number(process.env.PARACLIN_PORT) || 5199,
    proxy: {
      "/api": {
        target: process.env.PARACLIN_BACKEND || "http://127.0.0.1:8077",
        changeOrigin: true,
      },
    },
  },
});
