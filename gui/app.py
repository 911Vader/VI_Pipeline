"""
app.py -- VI_Pipeline desktop GUI  ·  Dark Lab Edition

Aesthetic: high-end genomics meets forensic intelligence.
Dark backgrounds, clinical precision, teal accents.
"""

from __future__ import annotations

import math
import queue
import sys
import tkinter as tk
import tkinter.font as tkfont
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext

import ttkbootstrap as tb
from tkinter import ttk

from gui.runner import start_pipeline_thread
from gui.validation import validate_inputs

# ── Palette ─────────────────────────────────────────────────────────────────

_BG      = "#0D1117"
_SURFACE = "#161B22"
_BORDER  = "#21262D"
_ACCENT  = "#00D4AA"
_ACCENT2 = "#1F6FEB"
_TEXT    = "#E6EDF3"
_MUTED   = "#8B949E"
_DANGER  = "#DA3633"
_SUCCESS = "#238636"
_WARNING = "#D29922"

_FONT_UI   = "Ubuntu"
_FONT_MONO = "Ubuntu Mono"
_SZ        = 11
_SZS       = 9

_PAD  = 16
_PADS = 8

_PIPELINE_ROOT = Path(__file__).parent.parent.resolve()
_LOGO_PATH     = Path(__file__).parent / "assets" / "logo.png"
_GOV_LOGO_PATH = Path(__file__).parent / "assets" / "dubai_gov.svg"

_TITLE = "VI_Pipeline  ·  Forensic Variant Interpretation  ·  Genome Center"
_STEPS = 7


# ── Helpers ──────────────────────────────────────────────────────────────────

def _hex_fade(hex_color: str, factor: float) -> str:
    """Darken a hex colour by factor (0–1)."""
    c = hex_color.lstrip("#")
    r, g, b = int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)
    r = int(r * factor); g = int(g * factor); b = int(b * factor)
    return f"#{r:02x}{g:02x}{b:02x}"


# ── Main window ──────────────────────────────────────────────────────────────

class App(tb.Window):

    def __init__(self) -> None:
        super().__init__(title=_TITLE, themename="darkly", minsize=(920, 720))
        self.configure(background=_BG)
        self._configure_fonts()
        self._apply_styles()

        self._queue: queue.Queue = queue.Queue()
        self._running: bool = False
        self._stopping: bool = False
        self._output_dir: str = ""
        self._proc_holder: list = [None]
        self._filter_expanded: bool = False
        self._filter_running: bool = False
        self._filter_queue: queue.Queue = queue.Queue()
        self._filter_proc: list = [None]
        self._run_start: float | None = None
        self._timer_id: str | None = None
        self._skip_vep_var: tk.BooleanVar | None = None
        self._second_stage_cmd: list | None = None
        self._pending_subset_values: dict | None = None
        self._pending_subset_out: str = ""
        self._pending_custom_filter_bash: str = ""   # custom filter run after vep_filter.sh
        self._stage_offset: int = 0   # added to ##STAGE## n when subset follows CLI
        self._last_run_values: dict | None = None    # saved for Resume
        self._last_run_out: str = ""                 # saved for Resume

        # stepper state
        self._step_current: int = 0
        self._step_total: int = _STEPS

        # pulse state
        self._pulse_id: str | None = None
        self._pulse_on: bool = True

        self._build_ui()
        self._set_status("Ready")

    # ── Fonts & styles ────────────────────────────────────────────────────────

    def _configure_fonts(self) -> None:
        try:
            tkfont.nametofont("TkDefaultFont").configure(family=_FONT_UI, size=_SZ)
            tkfont.nametofont("TkTextFont").configure(family=_FONT_UI, size=_SZ)
            tkfont.nametofont("TkFixedFont").configure(family=_FONT_MONO, size=_SZ)
        except tk.TclError:
            pass

    def _apply_styles(self) -> None:
        s = self.style
        s.configure(".", background=_BG, foreground=_TEXT, font=(_FONT_UI, _SZ))
        s.configure("TFrame",      background=_BG)
        s.configure("Surface.TFrame", background=_SURFACE)
        s.configure("TLabel",      background=_BG, foreground=_TEXT)
        s.configure("Surface.TLabel", background=_SURFACE, foreground=_TEXT)
        s.configure("Muted.TLabel",   background=_SURFACE, foreground=_MUTED,
                    font=(_FONT_UI, _SZS))
        s.configure("SectionHdr.TLabel", background=_SURFACE, foreground=_MUTED,
                    font=(_FONT_UI, _SZS, "bold"))
        s.configure("FieldLbl.TLabel",   background=_SURFACE, foreground=_TEXT,
                    font=(_FONT_UI, _SZ))
        s.configure("Accent.TLabel", background=_BG, foreground=_ACCENT,
                    font=(_FONT_UI, _SZ, "bold"))
        s.configure("Status.TLabel", background=_BG, foreground=_MUTED,
                    font=(_FONT_UI, _SZ, "italic"))
        s.configure("Elapsed.TLabel", background=_BG, foreground=_MUTED,
                    font=(_FONT_MONO, _SZ))
        s.configure("LogHdr.TLabel",  background=_SURFACE, foreground=_MUTED,
                    font=(_FONT_UI, _SZS, "bold"))
        s.configure("TEntry",
                    fieldbackground=_BG, foreground=_TEXT,
                    insertcolor=_ACCENT, font=(_FONT_UI, _SZ),
                    bordercolor=_BORDER, lightcolor=_BORDER, darkcolor=_BORDER,
                    padding=(6, 6))
        s.configure("TCheckbutton",
                    background=_SURFACE, foreground=_MUTED,
                    font=(_FONT_UI, _SZ))
        s.map("TCheckbutton",
              foreground=[("active", _TEXT)],
              background=[("active", _SURFACE)])
        s.configure("Teal.Horizontal.TProgressbar",
                    troughcolor=_BORDER, background=_ACCENT,
                    thickness=8, borderwidth=0, relief="flat")

    # ── Header ────────────────────────────────────────────────────────────────

    def _build_header(self, parent: tk.Frame) -> None:
        from PIL import Image, ImageTk
        import io

        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=2)
        parent.columnconfigure(2, weight=1)

        surf_rgb = (22, 27, 34)   # _SURFACE as RGB tuple

        def _load_png(path: Path, height: int):
            """Load PNG onto dark surface; removes near-white background pixels."""
            try:
                img = Image.open(path).convert("RGBA")
                pixels = img.load()
                for y in range(img.height):
                    for x in range(img.width):
                        r, g, b, a = pixels[x, y]
                        if r > 230 and g > 230 and b > 230:
                            pixels[x, y] = (r, g, b, 0)
                w = int(height * img.width / img.height)
                img = img.resize((w, height), Image.LANCZOS)
                bg = Image.new("RGBA", img.size, surf_rgb + (255,))
                bg.paste(img, mask=img.split()[3])
                return ImageTk.PhotoImage(bg.convert("RGB"))
            except Exception:
                return None

        def _load_svg_colorized(path: Path, height: int):
            """Load SVG composited on dark bg; dark text → white, emblem colors preserved."""
            try:
                import cairosvg
                svg_text = path.read_text(encoding="utf-8")
                # cls-4 is the near-black text (#19171c); make it white for dark bg
                svg_text = svg_text.replace("fill: #19171c", "fill: #ffffff")
                svg_text = svg_text.replace("fill:#19171c", "fill:#ffffff")
                png = cairosvg.svg2png(
                    bytestring=svg_text.encode(), output_width=800, output_height=800
                )
                img = Image.open(io.BytesIO(png)).convert("RGBA")
                bbox = img.getbbox()
                if bbox:
                    img = img.crop(bbox)
                w = int(height * img.width / img.height)
                img = img.resize((w, height), Image.LANCZOS)
                bg = Image.new("RGBA", img.size, surf_rgb + (255,))
                bg.paste(img, mask=img.split()[3])
                return ImageTk.PhotoImage(bg.convert("RGB"))
            except Exception:
                return None

        # ── left: Government of Dubai logo (white on dark)
        left = tk.Frame(parent, bg=_SURFACE)
        left.grid(row=0, column=0, sticky="nsew")
        gov_photo = _load_svg_colorized(_GOV_LOGO_PATH, 44)
        if gov_photo:
            lbl = tk.Label(left, image=gov_photo, bg=_SURFACE, padx=_PAD)
            lbl.image = gov_photo
            lbl.pack(side="left", pady=14)
        else:
            tk.Label(left, text="Government of Dubai",
                     bg=_SURFACE, fg="white",
                     font=(_FONT_UI, _SZS), padx=_PAD).pack(side="left", pady=14)

        # ── center: Genome Center title
        ctr = tk.Frame(parent, bg=_SURFACE)
        ctr.grid(row=0, column=1, sticky="nsew")
        tk.Label(ctr, text="G E N O M E   C E N T E R",
                 bg=_SURFACE, fg=_ACCENT,
                 font=(_FONT_UI, 15, "bold")).pack(pady=(16, 2))
        tk.Frame(ctr, bg=_ACCENT, height=2, width=160).pack()
        tk.Label(ctr, text="Forensic Variant Interpretation Pipeline",
                 bg=_SURFACE, fg=_MUTED,
                 font=(_FONT_UI, 9)).pack(pady=(4, 14))

        # ── right: Dubai Police logo
        right = tk.Frame(parent, bg=_SURFACE)
        right.grid(row=0, column=2, sticky="nsew")
        police_photo = _load_png(_LOGO_PATH, 68)
        if police_photo:
            lbl2 = tk.Label(right, image=police_photo, bg=_SURFACE, padx=_PAD)
            lbl2.image = police_photo
            lbl2.pack(side="right", pady=6)
        else:
            tk.Label(right, text="Dubai Police",
                     bg=_SURFACE, fg=_ACCENT,
                     font=(_FONT_UI, _SZS), padx=_PAD).pack(side="right", pady=14)

        # ── teal separator
        tk.Frame(parent, bg=_ACCENT, height=1).grid(
            row=1, column=0, columnspan=3, sticky="ew")

    # ── Layout skeleton ───────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(4, weight=1)

        # Header
        hdr = tk.Frame(self, bg=_SURFACE)
        hdr.grid(row=0, column=0, sticky="ew")
        self._build_header(hdr)

        # Form card
        self._build_form_card(row=1)

        # Controls
        ctrl = tk.Frame(self, bg=_BG)
        ctrl.grid(row=2, column=0, sticky="ew", padx=_PAD, pady=(_PADS, _PADS))
        ctrl.columnconfigure(0, weight=1)
        self._build_controls(ctrl)

        # Advanced filter card
        self._build_filter_card(row=3)

        # Log card
        self._build_log_card(row=4)

    # ── Form card ─────────────────────────────────────────────────────────────

    def _build_form_card(self, row: int) -> None:
        outer = tk.Frame(self, bg=_BORDER, padx=1, pady=1)
        outer.grid(row=row, column=0, sticky="ew", padx=_PAD, pady=(_PAD, _PADS))

        inner = tk.Frame(outer, bg=_SURFACE)
        inner.pack(fill="both", expand=True)
        inner.columnconfigure(2, weight=1)

        # 2 px teal left accent bar
        tk.Frame(inner, bg=_ACCENT, width=3).grid(
            row=0, column=0, rowspan=99, sticky="ns")

        content = tk.Frame(inner, bg=_SURFACE)
        content.grid(row=0, column=1, sticky="nsew", padx=(0, 0))
        content.columnconfigure(1, weight=1)
        inner.columnconfigure(1, weight=1)

        # Section header
        tk.Label(content, text="SAMPLE INFORMATION",
                 bg=_SURFACE, fg=_MUTED,
                 font=(_FONT_UI, _SZS, "bold"), anchor="w").grid(
            row=0, column=0, columnspan=4, sticky="w", padx=_PAD, pady=(_PAD, _PADS))

        tk.Frame(content, bg=_BORDER, height=1).grid(
            row=1, column=0, columnspan=4, sticky="ew", padx=_PAD)

        _ROWS = [
            ("input",     "Sample file",              "file",
             "VCF, FASTQ, or BAM — the sequencing file for this sample"),
            ("input2",    "Second reads  (optional)", "file",
             "R2 FASTQ only — leave blank for VCF / BAM input"),
            ("gene_list", "Gene list",                "file",
             "One gene name per line  (.txt)"),
            ("sample_id", "Sample ID",                None,
             "Short label used to name the output folder"),
            ("case_id",   "Case ID  (optional)",      None,
             "Investigation or case reference number"),
        ]

        self._vars: dict[str, tk.StringVar] = {}
        self._entry_frames: dict[str, tk.Frame] = {}

        for i, (key, label, browse_kind, hint) in enumerate(_ROWS):
            r = i + 2

            tk.Label(content, text=label, bg=_SURFACE, fg=_TEXT,
                     font=(_FONT_UI, _SZ), anchor="w", width=22).grid(
                row=r, column=0, sticky="w", padx=(_PAD, _PADS), pady=(8, 3))

            var = tk.StringVar()
            self._vars[key] = var

            # Entry with focus-highlight border frame
            ef = tk.Frame(content, bg=_BORDER, padx=1, pady=1)
            ef.grid(row=r, column=1, sticky="ew", pady=(8, 3))
            self._entry_frames[key] = ef

            entry = tk.Entry(
                ef, textvariable=var,
                bg=_BG, fg=_TEXT, insertbackground=_ACCENT,
                selectbackground=_ACCENT, selectforeground=_BG,
                relief="flat", font=(_FONT_UI, _SZ),
                highlightthickness=0, bd=4,
            )
            entry.pack(fill="x")

            entry.bind("<FocusIn>",  lambda e, f=ef: f.configure(bg=_ACCENT))
            entry.bind("<FocusOut>", lambda e, f=ef: f.configure(bg=_BORDER))
            self._attach_entry_context_menu(entry)

            if browse_kind == "file":
                browse_btn = tk.Button(
                    content, text="Browse",
                    command=lambda k=key: self._browse_file(k),
                    bg=_BORDER, fg=_MUTED, activebackground=_BORDER,
                    activeforeground=_ACCENT, relief="flat",
                    font=(_FONT_UI, _SZ), cursor="hand2",
                    padx=10, pady=5, bd=0,
                )
                browse_btn.grid(row=r, column=2, padx=(_PADS, 0), pady=(8, 3))
            else:
                browse_btn = None

            if key == "input2":
                self._r2_entry = entry
                self._r2_btn   = browse_btn

            tk.Label(content, text=hint, bg=_SURFACE, fg=_MUTED,
                     font=(_FONT_UI, _SZS), anchor="w").grid(
                row=r, column=3, sticky="w", padx=(_PADS, _PAD), pady=(8, 3))

        # Skip VEP
        self._skip_vep_var = tk.BooleanVar(value=False)
        skip_row = len(_ROWS) + 2
        skip_cb = tk.Checkbutton(
            content,
            text="Skip annotation  (only if file is already annotated with VEP)",
            variable=self._skip_vep_var,
            bg=_SURFACE, fg=_MUTED, selectcolor=_BG,
            activebackground=_SURFACE, activeforeground=_TEXT,
            font=(_FONT_UI, _SZS), cursor="hand2",
        )
        skip_cb.grid(row=skip_row, column=0, columnspan=4,
                     sticky="w", padx=_PAD, pady=(_PADS, 2))

        self._disable_small_model_var = tk.BooleanVar(value=False)
        dsm_cb = tk.Checkbutton(
            content,
            text="Deep model only  (FASTQ input: disable DeepVariant small model — forces full CNN on every site)",
            variable=self._disable_small_model_var,
            bg=_SURFACE, fg=_MUTED, selectcolor=_BG,
            activebackground=_SURFACE, activeforeground=_TEXT,
            font=(_FONT_UI, _SZS), cursor="hand2",
        )
        dsm_cb.grid(row=skip_row + 1, column=0, columnspan=4,
                    sticky="w", padx=_PAD, pady=(2, _PAD))

        tk.Frame(content, bg=_SURFACE, height=4).grid(
            row=skip_row + 2, column=0, columnspan=4)

        self._vars["input"].trace_add("write", self._on_input_changed)
        self._vars["input"].trace_add("write", self._on_form_changed)
        self._vars["sample_id"].trace_add("write", self._on_form_changed)

    # ── Controls ─────────────────────────────────────────────────────────────

    def _build_controls(self, parent: tk.Frame) -> None:
        parent.columnconfigure(0, weight=1)

        # ── Buttons
        btn_row = tk.Frame(parent, bg=_BG)
        btn_row.grid(row=0, column=0, sticky="ew")
        btn_row.columnconfigure(0, weight=1)

        self._run_btn = tk.Button(
            btn_row,
            text="▶   START ANALYSIS",
            command=self._on_run,
            bg=_ACCENT, fg=_BG, activebackground=_hex_fade(_ACCENT, 0.85),
            activeforeground=_BG,
            font=(_FONT_UI, 13, "bold"), relief="flat", bd=0,
            cursor="hand2", padx=20, pady=13,
        )
        self._run_btn.grid(row=0, column=0, sticky="ew", padx=(0, 6))

        self._stop_btn = tk.Button(
            btn_row,
            text="■   STOP",
            command=self._on_stop,
            state="disabled",
            bg=_DANGER, fg=_TEXT, activebackground=_hex_fade(_DANGER, 0.85),
            activeforeground=_TEXT,
            font=(_FONT_UI, 12, "bold"), relief="flat", bd=0,
            cursor="hand2", padx=20, pady=13,
        )
        self._stop_btn.grid(row=0, column=1, padx=(0, 6))

        self._resume_btn = tk.Button(
            btn_row,
            text="⟳   RESUME",
            command=self._on_resume,
            state="disabled",
            bg=_BORDER, fg=_MUTED, activebackground=_hex_fade(_BORDER, 0.85),
            activeforeground=_TEXT,
            font=(_FONT_UI, 12, "bold"), relief="flat", bd=0,
            cursor="hand2", padx=20, pady=13,
        )
        self._resume_btn.grid(row=0, column=2, padx=(0, 6))

        self._reset_btn = tk.Button(
            btn_row,
            text="↺   NEW RUN",
            command=self._on_new_run,
            bg=_BORDER, fg=_MUTED, activebackground=_hex_fade(_BORDER, 0.85),
            activeforeground=_TEXT,
            font=(_FONT_UI, 12, "bold"), relief="flat", bd=0,
            cursor="hand2", padx=20, pady=13,
        )
        self._reset_btn.grid(row=0, column=3)

        # ── Status row with pulsing dot
        status_row = tk.Frame(parent, bg=_BG)
        status_row.grid(row=1, column=0, sticky="ew", pady=(8, 4))
        status_row.columnconfigure(1, weight=1)

        self._pulse_canvas = tk.Canvas(
            status_row, width=12, height=12, bg=_BG,
            highlightthickness=0)
        self._pulse_canvas.grid(row=0, column=0, padx=(2, 6))
        self._pulse_dot = self._pulse_canvas.create_oval(
            2, 2, 10, 10, fill=_BG, outline=_BG)

        self._status_var = tk.StringVar()
        tk.Label(status_row, textvariable=self._status_var,
                 bg=_BG, fg=_MUTED, font=(_FONT_UI, _SZ, "italic"),
                 anchor="w").grid(row=0, column=1, sticky="w")

        self._elapsed_var = tk.StringVar()
        tk.Label(status_row, textvariable=self._elapsed_var,
                 bg=_BG, fg=_MUTED, font=(_FONT_MONO, _SZ),
                 anchor="e").grid(row=0, column=2, sticky="e")

        # ── 7-dot stepper
        self._stepper_canvas = tk.Canvas(
            parent, bg=_BG, height=28, highlightthickness=0)
        self._stepper_canvas.grid(row=2, column=0, sticky="ew", pady=(4, 4))
        self._stepper_canvas.bind("<Configure>", lambda e: self._draw_stepper())

        # ── Progress bar
        self._progress_var = tk.DoubleVar(value=0)
        self._progress_bar = ttk.Progressbar(
            parent, variable=self._progress_var,
            maximum=100, style="Teal.Horizontal.TProgressbar",
        )
        self._progress_bar.grid(row=3, column=0, sticky="ew", pady=(0, 4))

        # ── Stage label
        self._stage_var = tk.StringVar()
        tk.Label(parent, textvariable=self._stage_var,
                 bg=_BG, fg=_ACCENT, font=(_FONT_UI, _SZS, "bold"),
                 anchor="w").grid(row=4, column=0, sticky="w")

    # ── Stepper dots ──────────────────────────────────────────────────────────

    def _draw_stepper(self) -> None:
        c = self._stepper_canvas
        c.delete("all")
        W = c.winfo_width()
        if W < 10:
            return

        n = _STEPS
        r = 5            # dot radius
        cy = 14          # vertical centre
        margin = 24
        spacing = (W - 2 * margin) / max(n - 1, 1)

        for i in range(n):
            x = margin + i * spacing

            # connecting line
            if i < n - 1:
                nx = margin + (i + 1) * spacing
                if i < self._step_current - 1:
                    col = _ACCENT
                elif i == self._step_current - 1:
                    col = _hex_fade(_ACCENT, 0.5)
                else:
                    col = _BORDER
                c.create_line(x + r, cy, nx - r, cy, fill=col, width=2)

            # dot
            if i < self._step_current:
                # completed
                c.create_oval(x - r, cy - r, x + r, cy + r,
                              fill=_ACCENT, outline=_ACCENT)
                c.create_text(x, cy, text="✓", fill=_BG,
                              font=(_FONT_UI, 7, "bold"))
            elif i == self._step_current and self._running:
                # current — larger ring
                c.create_oval(x - r - 2, cy - r - 2, x + r + 2, cy + r + 2,
                              fill=_BG, outline=_ACCENT, width=2)
            else:
                # future
                c.create_oval(x - r, cy - r, x + r, cy + r,
                              fill=_BG, outline=_BORDER, width=1)

    # ── Pulsing dot ───────────────────────────────────────────────────────────

    def _pulse_start(self) -> None:
        self._pulse_on = True
        self._pulse_tick()

    def _pulse_stop(self) -> None:
        if self._pulse_id:
            self.after_cancel(self._pulse_id)
            self._pulse_id = None
        self._pulse_canvas.itemconfigure(self._pulse_dot, fill=_BG, outline=_BG)

    def _pulse_tick(self) -> None:
        if not self._running:
            self._pulse_stop()
            return
        color = _ACCENT if self._pulse_on else _hex_fade(_ACCENT, 0.25)
        self._pulse_canvas.itemconfigure(self._pulse_dot,
                                         fill=color, outline=color)
        self._pulse_on = not self._pulse_on
        self._pulse_id = self.after(600, self._pulse_tick)

    # ── Advanced subset filter ────────────────────────────────────────────────

    def _build_filter_card(self, row: int) -> None:
        outer = tk.Frame(self, bg=_BORDER, padx=1, pady=1)
        outer.grid(row=row, column=0, sticky="ew", padx=_PAD, pady=(0, _PADS))
        outer.columnconfigure(0, weight=1)

        inner = tk.Frame(outer, bg=_SURFACE)
        inner.pack(fill="both", expand=True)
        inner.columnconfigure(0, weight=1)

        # ── Toggle header
        hdr = tk.Frame(inner, bg=_SURFACE, cursor="hand2")
        hdr.pack(fill="x")
        tk.Frame(hdr, bg=_ACCENT2, width=3).pack(side="left", fill="y")
        self._filter_arrow = tk.Label(hdr, text="▶", bg=_SURFACE, fg=_ACCENT2,
                                      font=(_FONT_UI, _SZS), padx=6, cursor="hand2")
        self._filter_arrow.pack(side="left")
        tk.Label(hdr, text="ADVANCED SUBSET FILTER",
                 bg=_SURFACE, fg=_MUTED,
                 font=(_FONT_UI, _SZS, "bold"), pady=_PADS,
                 cursor="hand2").pack(side="left")

        for w in (hdr,) + tuple(hdr.winfo_children()):
            w.bind("<Button-1>", lambda e: self._toggle_filter())

        # ── Collapsible body
        self._filter_body = tk.Frame(inner, bg=_SURFACE)
        self._filter_body.columnconfigure(3, weight=1)

        # Separator
        tk.Frame(self._filter_body, bg=_BORDER, height=1).grid(
            row=0, column=0, columnspan=6, sticky="ew", padx=_PAD)

        # Match-mode row
        mf = tk.Frame(self._filter_body, bg=_SURFACE)
        mf.grid(row=1, column=0, columnspan=6, sticky="w", padx=_PAD, pady=(_PADS, 2))
        tk.Label(mf, text="Match:", bg=_SURFACE, fg=_MUTED,
                 font=(_FONT_UI, _SZS)).pack(side="left")
        self._filter_mode = tk.StringVar(value="ALL")
        for txt, val in [("ALL criteria  (AND)", "ALL"), ("ANY criterion  (OR)", "ANY")]:
            tk.Radiobutton(mf, text=txt, variable=self._filter_mode, value=val,
                           bg=_SURFACE, fg=_TEXT, selectcolor=_BG,
                           activebackground=_SURFACE, activeforeground=_TEXT,
                           font=(_FONT_UI, _SZS)).pack(side="left", padx=(10, 0))

        # Column headers
        for col, txt in enumerate(["", "Criterion", "Op", "Value"]):
            tk.Label(self._filter_body, text=txt, bg=_SURFACE, fg=_MUTED,
                     font=(_FONT_UI, _SZS, "bold"), anchor="w").grid(
                row=2, column=col, padx=(4, 2), pady=(2, 0), sticky="w")

        # Criteria definitions: (key, label, op, default, kind)
        # kind: "float" = numeric entry, "bool" = checkbox only
        self._fvars: dict = {}
        _CRIT = [
            ("cadd",     "CADD phred",    "≥", "20",    "float"),
            ("revel",    "REVEL",         "≥", "0.75",  "float"),
            ("spliceai", "SpliceAI max",  "≥", "0.5",   "float"),
            ("am",       "AlphaMissense", "≥", "0.564", "float"),
            ("af",       "gnomAD AF",     "≤", "0.001", "float"),
            ("loftee",   "LOFTEE  HC only", "", "",     "bool"),
            ("canon",    "Canonical only", "", "",      "bool"),
            ("genelist", "In gene list",  "", "",       "bool"),
        ]
        for i, (key, label, op, default, kind) in enumerate(_CRIT):
            r = 3 + i
            en = tk.BooleanVar(value=False)
            self._fvars[f"{key}_en"] = en
            tk.Checkbutton(self._filter_body, variable=en, bg=_SURFACE,
                           selectcolor=_BG, activebackground=_SURFACE,
                           cursor="hand2").grid(row=r, column=0, padx=(_PAD, 2),
                                                pady=1, sticky="w")
            tk.Label(self._filter_body, text=label, bg=_SURFACE, fg=_TEXT,
                     font=(_FONT_UI, _SZ), width=17, anchor="w").grid(
                row=r, column=1, padx=2, pady=1, sticky="w")
            tk.Label(self._filter_body, text=op, bg=_SURFACE, fg=_ACCENT,
                     font=(_FONT_MONO, _SZ, "bold"), width=2).grid(
                row=r, column=2, padx=2, pady=1)
            if kind == "float":
                vv = tk.StringVar(value=default)
                self._fvars[f"{key}_val"] = vv
                ef = tk.Frame(self._filter_body, bg=_BORDER, padx=1, pady=1)
                ef.grid(row=r, column=3, padx=(2, _PAD), pady=1, sticky="w")
                tk.Entry(ef, textvariable=vv, width=9,
                         bg=_BG, fg=_TEXT, insertbackground=_ACCENT,
                         relief="flat", font=(_FONT_MONO, _SZ),
                         bd=3, highlightthickness=0).pack()

        # ClinVar row
        cv_r = 3 + len(_CRIT)
        cv_en = tk.BooleanVar(value=False)
        self._fvars["clinvar_en"] = cv_en
        tk.Checkbutton(self._filter_body, variable=cv_en, bg=_SURFACE,
                       selectcolor=_BG, activebackground=_SURFACE,
                       cursor="hand2").grid(row=cv_r, column=0, padx=(_PAD, 2),
                                            pady=1, sticky="w")
        tk.Label(self._filter_body, text="ClinVar", bg=_SURFACE, fg=_TEXT,
                 font=(_FONT_UI, _SZ), width=17, anchor="w").grid(
            row=cv_r, column=1, padx=2, pady=1, sticky="w")
        cv_val = tk.StringVar(value="P / LP / Conflicting")
        self._fvars["clinvar_val"] = cv_val
        ttk.Combobox(self._filter_body, textvariable=cv_val, state="readonly",
                     width=26, font=(_FONT_UI, _SZ),
                     values=["P / LP / Conflicting",
                             "P / LP only",
                             "Conflicting only"]).grid(
            row=cv_r, column=2, columnspan=2, padx=(2, _PAD), pady=1, sticky="w")

        # Output name + apply button
        sep_r = cv_r + 1
        tk.Frame(self._filter_body, bg=_BORDER, height=1).grid(
            row=sep_r, column=0, columnspan=6, sticky="ew",
            padx=_PAD, pady=(_PADS, _PADS))

        out_r = sep_r + 1
        tk.Label(self._filter_body, text="Output filename",
                 bg=_SURFACE, fg=_MUTED,
                 font=(_FONT_UI, _SZS)).grid(
            row=out_r, column=0, columnspan=2, sticky="w", padx=_PAD)
        out_var = tk.StringVar(value="custom_filter")
        self._fvars["outname"] = out_var
        ef2 = tk.Frame(self._filter_body, bg=_BORDER, padx=1, pady=1)
        ef2.grid(row=out_r, column=2, columnspan=3, padx=(2, _PAD),
                 pady=2, sticky="ew")
        tk.Entry(ef2, textvariable=out_var, bg=_BG, fg=_TEXT,
                 insertbackground=_ACCENT, relief="flat",
                 font=(_FONT_UI, _SZ), bd=3, highlightthickness=0).pack(fill="x")

        btn_r = out_r + 1
        bf = tk.Frame(self._filter_body, bg=_SURFACE)
        bf.grid(row=btn_r, column=0, columnspan=6, sticky="w",
                padx=_PAD, pady=(_PADS, _PAD))
        self._filter_apply_btn = tk.Button(
            bf, text="⚡  Apply Filter",
            command=self._apply_custom_filter,
            bg=_ACCENT2, fg=_TEXT,
            activebackground=_hex_fade(_ACCENT2, 0.8),
            activeforeground=_TEXT,
            font=(_FONT_UI, _SZ, "bold"), relief="flat", bd=0,
            cursor="hand2", padx=14, pady=8,
        )
        self._filter_apply_btn.pack(side="left")
        self._filter_status_var = tk.StringVar(value="")
        tk.Label(bf, textvariable=self._filter_status_var,
                 bg=_SURFACE, fg=_MUTED,
                 font=(_FONT_UI, _SZS, "italic")).pack(side="left", padx=_PAD)

    def _toggle_filter(self) -> None:
        if self._filter_expanded:
            self._filter_body.pack_forget()
            self._filter_arrow.configure(text="▶")
            self._filter_expanded = False
        else:
            self._filter_body.pack(fill="both", expand=True)
            self._filter_arrow.configure(text="▼")
            self._filter_expanded = True

    def _build_custom_filter_bash(
        self,
        rare_vcf: str,
        subsets_dir: str,
    ) -> tuple[str, str] | tuple[None, str]:
        """
        Build the bcftools bash command for the custom subset filter.

        Returns (bash_cmd, filter_expr) on success, or (None, error_message) if
        no criteria are enabled or a required field (gene list) is missing.

        rare_vcf     — path to the input VCF (rare_annotated.1pc.vcf.gz)
        subsets_dir  — directory where the output .vcf.gz will be written
        """
        import tempfile as _tempfile

        out_name = self._fvars["outname"].get().strip() or "customized"
        if not out_name.endswith(".vcf.gz"):
            out_name += ".vcf.gz"
        out_vcf = str(Path(subsets_dir) / out_name)

        joiner = " && " if self._filter_mode.get() == "ALL" else " || "
        csq_fields = ["CANONICAL:String", "SYMBOL:String"]
        conditions: list[str] = []
        _tmp_gene: str = ""

        def _add(field: str, expr: str) -> None:
            if field and field not in csq_fields:
                csq_fields.append(field)
            conditions.append(expr)

        if self._fvars["cadd_en"].get():
            v = self._fvars["cadd_val"].get().strip()
            if v:
                _add("CADD_PHRED:Float", f"CADD_PHRED >= {v}")

        if self._fvars["revel_en"].get():
            v = self._fvars["revel_val"].get().strip()
            if v:
                _add("REVEL:Float", f"REVEL >= {v}")

        if self._fvars["spliceai_en"].get():
            v = self._fvars["spliceai_val"].get().strip()
            if v:
                _add("", (f"(SpliceAI_pred_DS_AG>={v} || SpliceAI_pred_DS_AL>={v}"
                          f" || SpliceAI_pred_DS_DG>={v} || SpliceAI_pred_DS_DL>={v})"))

        if self._fvars["am_en"].get():
            v = self._fvars["am_val"].get().strip()
            if v:
                _add("am_pathogenicity:Float", f"am_pathogenicity >= {v}")

        if self._fvars["af_en"].get():
            v = self._fvars["af_val"].get().strip()
            if v:
                _add("", f"MAX_AF <= {v}")

        if self._fvars["loftee_en"].get():
            _add("LoF:String", 'LoF="HC"')

        if self._fvars["canon_en"].get():
            _add("", "CANONICAL='YES'")

        if self._fvars["genelist_en"].get():
            gl = self._vars.get("gene_list", tk.StringVar()).get().strip()
            if not gl or not Path(gl).exists():
                return None, "Gene list path is not set or file not found."
            tf = _tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
            with open(gl) as gf:
                for ln in gf:
                    tf.write(ln.strip().upper() + "\n")
            tf.close()
            _tmp_gene = tf.name
            _add("", f"SYMBOL=@{_tmp_gene}")

        if self._fvars["clinvar_en"].get():
            cv = self._fvars["clinvar_val"].get()
            if cv == "P / LP / Conflicting":
                _add("ClinVar_CLNSIG:String", 'ClinVar_CLNSIG~"[Pp]athogenic"')
            elif cv == "P / LP only":
                _add("ClinVar_CLNSIG:String",
                     'ClinVar_CLNSIG~"[Pp]athogenic" && ClinVar_CLNSIG!~"[Cc]onflicting"')
            elif cv == "Conflicting only":
                _add("ClinVar_CLNSIG:String", 'ClinVar_CLNSIG~"[Cc]onflicting"')

        if not conditions:
            return None, ""

        filter_expr = joiner.join(f"({c})" for c in conditions)
        cleanup = f"rm -f '{_tmp_gene}'" if _tmp_gene else "true"
        bash_cmd = (
            f"bcftools +split-vep '{rare_vcf}'"
            f" -c {','.join(csq_fields)}"
            f" -i '{filter_expr}'"
            f" -Oz -o '{out_vcf}'"
            f" && tabix '{out_vcf}'"
            f" && echo \"Customized filter done — $(bcftools view -H '{out_vcf}' | wc -l)"
            f" variants → subsets/{out_name}\""
            f" ; {cleanup}"
        )
        return bash_cmd, filter_expr

    def _apply_custom_filter(self) -> None:
        """Standalone Apply — runs the custom filter immediately on an existing VCF."""
        rare_vcf = self._locate_rare_vcf()
        if rare_vcf is None:
            messagebox.showerror(
                "File not found",
                "Cannot locate a rare-variant VCF to filter.\n\n"
                "Either:\n"
                "  • Run the full pipeline first, or\n"
                "  • Set the Sample file field to a VCF / rare_annotated.1pc.vcf.gz "
                "from a previous run.",
            )
            return
        if self._filter_running:
            return

        # Work out where to write: same folder as the source VCF's subsets/
        rare_path = Path(rare_vcf)
        subsets_dir = rare_path.parent / "subsets"
        subsets_dir.mkdir(parents=True, exist_ok=True)

        bash_cmd, filter_expr = self._build_custom_filter_bash(rare_vcf, str(subsets_dir))
        if bash_cmd is None:
            if filter_expr:
                messagebox.showwarning("Gene list", filter_expr)
            else:
                messagebox.showwarning("No criteria", "Enable at least one filter criterion.")
            return

        self._filter_running = True
        self._filter_apply_btn.configure(state="disabled")
        self._filter_status_var.set("Running…")
        self._log_line("", tag=None)
        self._log_line("─" * 60, tag=None)
        self._log_line(f"[Filter]  {filter_expr}", tag="cmd")
        self._log_line(f"[Input]   {rare_vcf}", tag="cmd")
        self._log_line(f"[Output]  {subsets_dir}", tag="cmd")
        self._log_line("", tag=None)

        self._filter_queue = queue.Queue()
        self._filter_proc  = [None]
        start_pipeline_thread(["bash", "-c", bash_cmd],
                              self._filter_queue, self._filter_proc)
        self.after(100, self._poll_filter)

    def _poll_filter(self) -> None:
        try:
            while True:
                kind, val = self._filter_queue.get_nowait()
                if kind == "line":
                    self._log_line(val, tag=None)
                elif kind in ("exit", "error"):
                    self._filter_running = False
                    self._filter_apply_btn.configure(state="normal")
                    if kind == "exit" and val == 0:
                        self._filter_status_var.set("Done ✓")
                    else:
                        msg = val if kind == "error" else f"exit code {val}"
                        self._filter_status_var.set(f"Error — {msg}")
                        self._log_line(f"Filter error: {msg}", tag="error")
                    return
        except queue.Empty:
            pass
        self.after(100, self._poll_filter)

    # ── Log card ──────────────────────────────────────────────────────────────

    def _copy_log_selection(self) -> None:
        """Copy selected log text, or the entire log if nothing is selected."""
        try:
            text = self._log.get(tk.SEL_FIRST, tk.SEL_LAST)
        except tk.TclError:
            text = self._log.get("1.0", tk.END)
        self.clipboard_clear()
        self.clipboard_append(text)

    def _copy_log_all(self) -> None:
        """Copy the entire log to the clipboard."""
        text = self._log.get("1.0", tk.END)
        self.clipboard_clear()
        self.clipboard_append(text)

    def _build_log_card(self, row: int) -> None:
        outer = tk.Frame(self, bg=_BORDER, padx=1, pady=1)
        outer.grid(row=row, column=0, sticky="nsew", padx=_PAD, pady=(_PADS, _PAD))
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(1, weight=1)

        hdr = tk.Frame(outer, bg=_SURFACE)
        hdr.grid(row=0, column=0, sticky="ew")
        tk.Label(hdr, text="ANALYSIS LOG",
                 bg=_SURFACE, fg=_MUTED,
                 font=(_FONT_UI, _SZS, "bold"), anchor="w").pack(
            side="left", padx=_PAD, pady=(_PADS, _PADS))

        # Permanent Copy All button — always visible, no hover needed
        tk.Button(
            hdr, text="Copy All",
            command=self._copy_log_all,
            bg=_BORDER, fg=_MUTED, activebackground=_SURFACE,
            activeforeground=_ACCENT, relief="flat", bd=0,
            font=(_FONT_UI, _SZS), cursor="hand2", padx=8, pady=4,
        ).pack(side="right", padx=(_PADS, _PAD), pady=4)

        tk.Button(
            hdr, text="Select All",
            command=lambda: (self._log.focus_set(),
                             self._log.tag_add("sel", "1.0", tk.END)),
            bg=_BORDER, fg=_MUTED, activebackground=_SURFACE,
            activeforeground=_ACCENT, relief="flat", bd=0,
            font=(_FONT_UI, _SZS), cursor="hand2", padx=8, pady=4,
        ).pack(side="right", pady=4)

        self._log = scrolledtext.ScrolledText(
            outer,
            wrap=tk.WORD,
            font=(_FONT_MONO, _SZ),
            bg=_BG, fg=_MUTED,
            insertbackground=_ACCENT,
            selectbackground=_ACCENT2,
            selectforeground=_TEXT,
            height=14,
            relief="flat", bd=0,
            padx=12, pady=8,
        )
        self._log.grid(row=1, column=0, sticky="nsew")

        self._log.tag_config("error", foreground=_DANGER)
        self._log.tag_config("ok",    foreground=_SUCCESS)
        self._log.tag_config("info",  foreground=_ACCENT)
        self._log.tag_config("warn",  foreground=_WARNING)
        self._log.tag_config("cmd",   foreground=_ACCENT2)

        _NAV = {"Up", "Down", "Left", "Right", "Home", "End",
                "Prior", "Next", "Control_L", "Control_R",
                "Shift_L", "Shift_R", "Alt_L", "Alt_R", "caps_lock", "Tab"}

        def _block(event: tk.Event) -> str | None:
            if event.state & 0x4:
                return None
            if event.keysym in _NAV:
                return None
            return "break"

        self._log.bind("<Key>",       _block)
        self._log.bind("<<Paste>>",   lambda e: "break")
        self._log.bind("<Control-c>", lambda e: (self._copy_log_selection(), "break"))
        self._log.bind("<Control-C>", lambda e: (self._copy_log_selection(), "break"))
        self._log.bind("<Control-a>", lambda e: (self._log.tag_add("sel", "1.0", tk.END), "break"))
        self._log.bind("<Control-A>", lambda e: (self._log.tag_add("sel", "1.0", tk.END), "break"))
        self._attach_log_context_menu()

    # ── Context menus ─────────────────────────────────────────────────────────

    def _attach_entry_context_menu(self, entry: tk.Entry) -> None:
        menu = tk.Menu(self, tearoff=0, bg=_SURFACE, fg=_TEXT,
                       activebackground=_ACCENT, activeforeground=_BG,
                       font=(_FONT_UI, _SZ), bd=1, relief="flat")
        menu.add_command(label="Cut",
                         command=lambda: (entry.focus_set(),
                                          entry.event_generate("<<Cut>>")))
        menu.add_command(label="Copy",
                         command=lambda: (entry.focus_set(),
                                          entry.event_generate("<<Copy>>")))
        menu.add_command(label="Paste",
                         command=lambda: (entry.focus_set(),
                                          entry.event_generate("<<Paste>>")))
        menu.add_separator()
        menu.add_command(label="Select All",
                         command=lambda: (entry.focus_set(),
                                          entry.selection_range(0, tk.END)))
        entry.bind("<Button-3>",
                   lambda e: menu.tk_popup(e.x_root, e.y_root))

    def _attach_log_context_menu(self) -> None:
        menu = tk.Menu(self._log, tearoff=0, bg=_SURFACE, fg=_TEXT,
                       activebackground=_ACCENT, activeforeground=_BG,
                       font=(_FONT_UI, _SZ), bd=1, relief="flat")
        menu.add_command(label="Copy",      command=self._copy_log_selection)
        menu.add_command(label="Copy All",  command=self._copy_log_all)
        menu.add_command(label="Select All",
                         command=lambda: self._log.tag_add("sel", "1.0", tk.END))

        def _show_menu(event: tk.Event) -> None:
            # Preserve existing selection before the popup opens
            try:
                sel_start = self._log.index(tk.SEL_FIRST)
                sel_end   = self._log.index(tk.SEL_LAST)
            except tk.TclError:
                sel_start = sel_end = None
            menu.tk_popup(event.x_root, event.y_root)
            if sel_start and sel_end:
                self._log.tag_add("sel", sel_start, sel_end)

        self._log.bind("<Button-3>", _show_menu)

    # ── Input detection ───────────────────────────────────────────────────────

    _FASTQ_EXTS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")

    def _on_input_changed(self, *_) -> None:
        path = self._vars["input"].get().strip().lower()
        is_fastq = any(path.endswith(e) for e in self._FASTQ_EXTS)
        state = "normal" if is_fastq else "disabled"
        self._r2_entry.configure(state=state)
        if self._r2_btn:
            self._r2_btn.configure(state=state)
        if not is_fastq:
            self._vars["input2"].set("")

    def _on_form_changed(self, *_) -> None:
        """Enable Resume whenever the form has enough content and no run is active."""
        if self._running:
            return
        has_input = bool(self._vars["input"].get().strip())
        has_id    = bool(self._vars["sample_id"].get().strip())
        state = "normal" if (has_input and has_id) else "disabled"
        self._resume_btn.configure(state=state)

    # ── File browsing ─────────────────────────────────────────────────────────

    _BROWSE_TYPES: dict = {
        "input": [
            ("VCF / FASTQ / BAM",
             "*.vcf *.vcf.gz *.bcf *.fastq *.fastq.gz *.fq *.fq.gz *.bam *.cram"),
            ("VCF files",        "*.vcf *.vcf.gz *.bcf"),
            ("FASTQ files",      "*.fastq *.fastq.gz *.fq *.fq.gz"),
            ("BAM / CRAM files", "*.bam *.cram"),
            ("All files", "*"),
        ],
        "input2": [("FASTQ files", "*.fastq *.fastq.gz *.fq *.fq.gz"), ("All", "*")],
        "gene_list": [("Text files", "*.txt"), ("All files", "*")],
    }

    def _browse_file(self, key: str) -> None:
        current = self._vars[key].get().strip()
        init = (str(Path(current).parent) if current and Path(current).exists()
                else str(_PIPELINE_ROOT / "config") if key == "gene_list"
                else str(_PIPELINE_ROOT))
        path = filedialog.askopenfilename(
            title=f"Select {key.replace('_', ' ')}",
            filetypes=self._BROWSE_TYPES.get(key, [("All files", "*")]),
            initialdir=init,
        )
        if path:
            if path.startswith("recent:///"):
                path = self._resolve_recent_uri(path)
                if not path:
                    messagebox.showwarning(
                        "Recent files not supported",
                        "Please navigate to the file directly in the browser "
                        "(do not select from the Recent section in the left sidebar).",
                    )
                    return
            self._vars[key].set(path)

    def _resolve_recent_uri(self, uri: str) -> str:
        """Try to resolve a GTK recent:/// URI to an actual filesystem path."""
        import subprocess
        try:
            result = subprocess.run(
                ["gio", "info", "--attributes=standard::target-uri,standard::symlink-target", uri],
                capture_output=True, text=True, timeout=3,
            )
            for line in result.stdout.splitlines():
                for prefix in ("standard::target-uri:", "uri:"):
                    if line.strip().lower().startswith(prefix):
                        raw = line.split(":", 1)[1].strip()
                        if raw.startswith("file://"):
                            from urllib.parse import unquote
                            return unquote(raw[7:])
        except Exception:
            pass
        return ""

    # ── Timer ─────────────────────────────────────────────────────────────────

    def _start_timer(self) -> None:
        import time
        self._run_start = time.monotonic()
        self._tick_timer()

    def _stop_timer(self) -> None:
        if self._timer_id:
            self.after_cancel(self._timer_id)
            self._timer_id = None

    def _tick_timer(self) -> None:
        import time
        if not self._running or self._run_start is None:
            return
        e = int(time.monotonic() - self._run_start)
        h, rem = divmod(e, 3600)
        m, s = divmod(rem, 60)
        self._elapsed_var.set(f"  {h:02d}:{m:02d}:{s:02d}")
        self._timer_id = self.after(1000, self._tick_timer)

    # ── Stage progress ────────────────────────────────────────────────────────

    def _on_stage_line(self, line: str) -> None:
        try:
            rest = line[len("##STAGE##"):].strip()
            fraction, label = rest.split(" ", 1)
            n_str, total_str = fraction.split("/")
            n, total = int(n_str), int(total_str)
            # Offset applied when subset stage follows main pipeline
            n_abs     = self._stage_offset + n
            total_abs = self._stage_offset + total
            self._progress_var.set(n_abs / total_abs * 100)
            self._step_current = n_abs
            self._step_total   = total_abs
            self._stage_var.set(f"Step {n_abs} of {total_abs}  —  {label.strip()}")
            self._draw_stepper()
        except Exception:
            pass

    # ── Run ───────────────────────────────────────────────────────────────────

    def _start_run(self, cmd: list, output_dir: str, status_msg: str) -> None:
        """Fresh run — clears log, resets all progress state."""
        self._output_dir = output_dir
        self._run_btn.configure(state="disabled")
        self._stop_btn.configure(state="normal")
        self._running = True
        self._stopping = False
        self._stage_offset = 0
        self._proc_holder = [None]
        self._clear_log()
        self._progress_var.set(0)
        self._step_current = 0
        self._stage_var.set("")
        self._elapsed_var.set("")
        self._set_status(status_msg, running=True)
        self._draw_stepper()
        self._start_timer()
        self._pulse_start()
        self._log_line(f"$ {' '.join(cmd)}", tag="cmd")
        self._log_line("", tag=None)
        start_pipeline_thread(cmd, self._queue, self._proc_holder)
        self.after(100, self._poll)

    def _continue_run(self, cmd: list, status_msg: str) -> None:
        """Chain a second stage — appends to the existing log and continues progress."""
        self._run_btn.configure(state="disabled")
        self._stop_btn.configure(state="normal")
        self._running = True
        self._stopping = False
        # Offset so subset ##STAGE## 1/6 … 6/6 appears as step N+1 … N+6
        self._stage_offset = self._step_current
        self._proc_holder = [None]
        self._set_status(status_msg, running=True)
        self._pulse_start()
        self._log_line("", tag=None)
        self._log_line("─" * 60, tag=None)
        self._log_line(f"$ {' '.join(cmd)}", tag="cmd")
        self._log_line("", tag=None)
        start_pipeline_thread(cmd, self._queue, self._proc_holder)
        self.after(100, self._poll)

    def _on_run(self) -> None:
        values = {k: v.get() for k, v in self._vars.items()}
        ok, msg = validate_inputs(values)
        if not ok:
            messagebox.showerror("Cannot start analysis", msg, parent=self)
            return
        sample_id = values["sample_id"].strip()
        out = str(_PIPELINE_ROOT / "runs" / sample_id)
        self._last_run_values = values
        self._last_run_out    = out
        self._resume_btn.configure(state="disabled")
        # Store values so we can build the subset command after the first stage
        # finishes — we need to know which annotated VCF was produced first.
        self._pending_subset_values = values
        self._pending_subset_out = out
        self._second_stage_cmd = None   # built in _on_finished

        # Capture custom filter criteria now so they run automatically after
        # vep_filter.sh completes.  The rare VCF path is fixed: vep_filter.sh
        # always writes it to <out>/vep_filter/rare_annotated.1pc.vcf.gz.
        rare_vcf_path = str(Path(out) / "vep_filter" / "rare_annotated.1pc.vcf.gz")
        subsets_path  = str(Path(out) / "vep_filter" / "subsets")
        bash_cmd, _   = self._build_custom_filter_bash(rare_vcf_path, subsets_path)
        self._pending_custom_filter_bash = bash_cmd or ""

        self._start_run(self._build_command(values, out), out,
                        "Running — this may take several hours for a whole-genome sample")

    def _on_resume(self) -> None:
        if self._running:
            return
        # Use last-run values from this session, or fall back to current form state
        if self._last_run_values:
            values = self._last_run_values
            out    = self._last_run_out
        else:
            values = {k: v.get() for k, v in self._vars.items()}
            ok, msg = validate_inputs(values)
            if not ok:
                messagebox.showerror("Cannot resume", msg, parent=self)
                return
            sample_id = values["sample_id"].strip()
            out = str(_PIPELINE_ROOT / "runs" / sample_id)
        self._last_run_values = values
        self._last_run_out    = out
        self._resume_btn.configure(state="disabled")
        self._pending_subset_values = values
        self._pending_subset_out    = out
        self._second_stage_cmd      = None
        rare_vcf_path = str(Path(out) / "vep_filter" / "rare_annotated.1pc.vcf.gz")
        subsets_path  = str(Path(out) / "vep_filter" / "subsets")
        bash_cmd, _   = self._build_custom_filter_bash(rare_vcf_path, subsets_path)
        self._pending_custom_filter_bash = bash_cmd or ""
        self._start_run(self._build_command(values, out), out,
                        "Resuming — skipping completed steps...")

    def _on_stop(self) -> None:
        if not self._running or self._stopping:
            return
        proc = self._proc_holder[0]
        if proc is None:
            self._stopping = True
            self._set_status("Stopping...")
            return
        if not messagebox.askyesno(
            "Stop Analysis",
            "Stop the running analysis?\n\nPartial output files remain in the output folder.",
            parent=self,
        ):
            return
        self._stopping = True
        self._stop_btn.configure(state="disabled")
        self._set_status("Stopping — killing all processes...")
        self._log_line("", tag=None)
        self._log_line("Stop requested — killing pipeline and all child processes...",
                       tag="error")
        try:
            import os, signal
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

    def _on_new_run(self) -> None:
        """Reset the GUI to a clean state so the user can start a new analysis."""
        if self._running:
            messagebox.showwarning(
                "Analysis running",
                "Stop the current analysis before starting a new run.",
                parent=self,
            )
            return
        if not messagebox.askyesno(
            "New Run",
            "Clear all fields and the log to start a fresh analysis?",
            parent=self,
        ):
            return

        # Clear all input fields
        for var in self._vars.values():
            var.set("")

        # Reset state
        self._output_dir = ""
        self._running = False
        self._stopping = False
        self._filter_running = False
        self._second_stage_cmd = None
        self._pending_subset_values = None
        self._pending_subset_out = ""
        self._pending_custom_filter_bash = ""
        self._stage_offset = 0
        self._run_start = None
        self._last_run_values = None
        self._last_run_out = ""
        self._skip_vep_var.set(False)
        self._disable_small_model_var.set(False)
        self._resume_btn.configure(state="disabled")

        # Reset UI
        self._clear_log()
        self._progress_var.set(0)
        self._step_current = 0
        self._step_total = 0
        self._stage_var.set("")
        self._elapsed_var.set("")
        self._status_var.set("")
        self._stop_timer()
        self._pulse_stop()
        self._draw_stepper()
        self._run_btn.configure(state="normal")
        self._stop_btn.configure(state="disabled")
        if hasattr(self, "_filter_status_var"):
            self._filter_status_var.set("")

    def _build_command(self, values: dict, out: str) -> list:
        cli = str(Path(__file__).parent.parent / "cli.py")
        cmd = [sys.executable, cli,
               "--input",     values["input"].strip(),
               "--gene-list", values["gene_list"].strip(),
               "--output",    out,
               "--sample-id", values["sample_id"].strip()]
        if values.get("case_id", "").strip():
            cmd += ["--case-id", values["case_id"].strip()]
        if values.get("input2", "").strip():
            cmd += ["--input2", values["input2"].strip()]
        if self._disable_small_model_var and self._disable_small_model_var.get():
            cmd.append("--disable-small-model")
        return cmd

    def _find_annotated_vcf(self, out_dir: str) -> str | None:
        """Return the annotated VCF produced by the CLI, or None if not found."""
        p = Path(out_dir)
        for pattern in ["*.norm.vep.vcf", "*.vep.vcf", "annotated.vcf.gz", "*.vep.vcf.gz"]:
            matches = sorted(p.glob(pattern))
            if matches:
                return str(matches[0])
        return None

    def _locate_rare_vcf(self) -> str | None:
        """
        Find the rare-variant VCF to use as input for the Advanced Subset Filter.

        Search order:
          1. <output_dir>/rare_annotated.1pc.vcf.gz  (result of a fresh pipeline run)
          2. <parent of input file>/rare_annotated.1pc.vcf.gz  (previous run's folder)
          3. The input file itself, if it is a .vcf.gz or .vcf  (user pointed directly
             at the rare VCF or any annotated VCF)
        """
        # 1. Current output dir (set after a pipeline run in this session)
        if self._output_dir:
            candidate = Path(self._output_dir) / "rare_annotated.1pc.vcf.gz"
            if candidate.exists():
                return str(candidate)

        # 2. Parent directory of whatever the user has in the Sample file field
        input_path = self._vars["input"].get().strip()
        if input_path:
            parent = Path(input_path).parent
            candidate = parent / "rare_annotated.1pc.vcf.gz"
            if candidate.exists():
                return str(candidate)

            # 3. The input file itself if it's a VCF
            p = Path(input_path)
            if p.exists() and (input_path.endswith(".vcf.gz") or input_path.endswith(".vcf")):
                return str(p)

        return None

    def _build_vep_filter_command(self, values: dict, out: str,
                                   annotated_vcf: str | None = None) -> list:
        script = str(Path(__file__).parent.parent / "vep_filter.sh")
        subset_out = str(Path(out) / "vep_filter")
        # Use the annotated VCF if available — avoids re-running VEP
        input_vcf = annotated_vcf if annotated_vcf else values["input"].strip()
        cmd = ["bash", script,
               "--input",     input_vcf,
               "--out-dir",   subset_out,
               "--gene-list", values["gene_list"].strip()]
        # Always skip VEP when we have an already-annotated VCF
        if annotated_vcf or (self._skip_vep_var and self._skip_vep_var.get()):
            cmd.append("--skip-vep")
        return cmd

    # ── Polling ───────────────────────────────────────────────────────────────

    def _poll(self) -> None:
        try:
            while True:
                kind, payload = self._queue.get_nowait()
                if kind == "line":
                    if payload.startswith("##STAGE##"):
                        self._on_stage_line(payload)
                        self._log_line(f"  ▶  {self._stage_var.get()}", tag="info")
                    else:
                        tag = None
                        pl = payload.lower()
                        if any(w in pl for w in ("error", "failed", "traceback", "exception")):
                            tag = "error"
                        elif any(w in pl for w in ("warning", "warn")):
                            tag = "warn"
                        self._log_line(payload, tag=tag)
                elif kind == "exit":
                    self._on_finished(payload)
                    return
                elif kind == "error":
                    self._log_line(f"[GUI] {payload}", tag="error")
        except queue.Empty:
            pass
        if self._running:
            self.after(100, self._poll)

    def _on_finished(self, returncode: int) -> None:
        self._running = False
        self._stop_timer()
        self._pulse_stop()
        self._log_line("", tag=None)

        if self._stopping:
            self._second_stage_cmd = None
            self._pending_custom_filter_bash = ""
            self._run_btn.configure(state="normal")
            self._stop_btn.configure(state="disabled")
            if self._last_run_values:
                self._resume_btn.configure(state="normal")
            self._log_line("Analysis stopped by user.", tag="error")
            self._set_status("Stopped.  Click ⟳ RESUME to continue from where it stopped.")
            self._draw_stepper()
            return

        if returncode != 0:
            self._second_stage_cmd = None
            self._pending_custom_filter_bash = ""
            self._run_btn.configure(state="normal")
            self._stop_btn.configure(state="disabled")
            if self._last_run_values:
                self._resume_btn.configure(state="normal")
            self._log_line(
                f"Something went wrong (exit code {returncode}). Check the log above.",
                tag="error")
            self._set_status(
                "Failed — see the log for details.  Click ⟳ RESUME to retry from the last completed step.")
            self._draw_stepper()
            return

        # Build and chain the subset stage without resetting the log or progress.
        if hasattr(self, "_pending_subset_values") and self._pending_subset_values:
            values = self._pending_subset_values
            out    = self._pending_subset_out
            self._pending_subset_values = None
            annotated = self._find_annotated_vcf(self._output_dir)
            if annotated:
                self._log_line(f"  Annotated VCF: {annotated}", tag="info")
            else:
                self._log_line(
                    "  Warning: annotated VCF not found — subset step will re-annotate.",
                    tag="warn")
            cmd = self._build_vep_filter_command(values, out, annotated_vcf=annotated)
            self._log_line("Analysis complete — producing variant subsets...", tag="ok")
            self._continue_run(cmd, "Producing variant subsets...")
            return

        # ── Stage 3: customized subset filter (runs after vep_filter.sh if configured)
        if self._pending_custom_filter_bash:
            bash_cmd = self._pending_custom_filter_bash
            self._pending_custom_filter_bash = ""
            self._log_line("", tag=None)
            self._log_line("Running customized subset filter...", tag="ok")
            self._log_line(f"$ bash -c {bash_cmd[:120]}{'…' if len(bash_cmd) > 120 else ''}",
                           tag="cmd")
            self._continue_run(["bash", "-c", bash_cmd], "Applying customized subset filter…")
            return

        # ── All done
        self._run_btn.configure(state="normal")
        self._stop_btn.configure(state="disabled")
        self._progress_var.set(100)
        self._step_current = _STEPS
        self._stage_var.set("")
        self._draw_stepper()
        self._log_line("✓  Analysis complete.", tag="ok")
        self._set_status("✓  Analysis complete", done=True)
        self._elapsed_var.set(
            self._elapsed_var.get())   # keep elapsed shown

    # ── Log helpers ───────────────────────────────────────────────────────────

    def _log_line(self, text: str, *, tag: str | None = None) -> None:
        if tag:
            self._log.insert(tk.END, text + "\n", tag)
        else:
            self._log.insert(tk.END, text + "\n")
        self._log.see(tk.END)

    def _clear_log(self) -> None:
        self._log.delete("1.0", tk.END)

    def _set_status(self, text: str, *, running: bool = False,
                    done: bool = False) -> None:
        self._status_var.set(f"  {text}")
