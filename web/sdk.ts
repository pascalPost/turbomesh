/**
 * Lightweight TypeScript wrapper for the turbomesh WebAssembly module.
 * Exposes helpers to load the wasm, initialize the mesh, and read block point data.
 */

export type BlockSize = { i: number; j: number };

export type BlockPoints = {
  size: BlockSize;
  values: Float64Array; // view into wasm memory; copy if you need to hold it after calling freeMesh
};

export interface TurbomeshWasmExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  run(): void;
  freeMesh(): void;
  blocksCount(): number;
  blockSizeI(blockIdx: number): number;
  blockSizeJ(blockIdx: number): number;
  blockSize(blockIdx: number, outPtr: number): void; // legacy pointer-based helper
  blockPointsPtr(blockIdx: number): number;
  blockPointsLen(blockIdx: number): number;
}

export type LoaderOptions = {
  wasmUrl: string;
  /** Additional imports to pass through; `env.console_log` is provided if you do not override it. */
  imports?: WebAssembly.Imports;
  /** Optional log handler for `env.console_log` output. */
  onLog?: (message: string) => void;
  /** RequestInit forwarded to fetch; useful for credentials/custom headers. */
  fetchOptions?: RequestInit;
  /** Call `run()` immediately after instantiation. Defaults to true. */
  autoInit?: boolean;
};

const textDecoder = new TextDecoder("utf-8");

export class TurboMeshSDK {
  private constructor(private readonly exports: TurbomeshWasmExports) {}

  static async load(options: LoaderOptions): Promise<TurboMeshSDK> {
    const { wasmUrl, fetchOptions, autoInit = true, onLog } = options;
    let memory: WebAssembly.Memory | undefined;

    const userImports = options.imports ?? {};
    const userEnv = (userImports as { env?: Record<string, unknown> }).env ?? {};
    const logSink = onLog ?? ((message: string) => console.log(message));

    const defaultEnv = {
      console_log: (ptr: number, len: number) => {
        if (!memory) return;
        const bytes = new Uint8Array(memory.buffer, ptr, len);
        logSink(textDecoder.decode(bytes).trimEnd());
      },
    };

    const imports: WebAssembly.Imports = {
      ...userImports,
      env: { ...defaultEnv, ...userEnv },
    };

    const response = await fetch(wasmUrl, fetchOptions);
    let instance: WebAssembly.Instance;
    if ("instantiateStreaming" in WebAssembly && WebAssembly.instantiateStreaming) {
      try {
        const result = await WebAssembly.instantiateStreaming(response.clone(), imports);
        instance = result.instance;
      } catch {
        const buffer = await response.arrayBuffer();
        const result = await WebAssembly.instantiate(buffer, imports);
        instance = result.instance;
      }
    } else {
      const buffer = await response.arrayBuffer();
      const result = await WebAssembly.instantiate(buffer, imports);
      instance = result.instance;
    }

    const exports = instance.exports as unknown as TurbomeshWasmExports;
    memory = exports.memory;

    const sdk = new TurboMeshSDK(exports);
    if (autoInit) sdk.run();
    return sdk;
  }

  run(): void {
    this.exports.run();
  }

  free(): void {
    this.exports.freeMesh();
  }

  blocksCount(): number {
    return this.exports.blocksCount();
  }

  blockSize(blockIdx: number): BlockSize {
    return {
      i: this.exports.blockSizeI(blockIdx),
      j: this.exports.blockSizeJ(blockIdx),
    };
  }

  blockPointsView(blockIdx: number): BlockPoints {
    const ptr = this.exports.blockPointsPtr(blockIdx);
    const len = this.exports.blockPointsLen(blockIdx);
    if (ptr === 0 || len === 0) {
      throw new Error(`No point data for block ${blockIdx}`);
    }

    return {
      size: this.blockSize(blockIdx),
      values: new Float64Array(this.exports.memory.buffer, ptr, len),
    };
  }

  blockPointsCopy(blockIdx: number): BlockPoints {
    const view = this.blockPointsView(blockIdx);
    return { size: view.size, values: new Float64Array(view.values) };
  }
}
