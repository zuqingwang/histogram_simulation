# 执行 notebook 并把结果写出（nbconvert CLI 不可用，直接用 nbclient）
import sys, time
import nbformat
from nbclient import NotebookClient

src = sys.argv[1]
dst = sys.argv[2] if len(sys.argv) > 2 else src.replace(".ipynb", "_out.ipynb")
timeout = int(sys.argv[3]) if len(sys.argv) > 3 else 7200

nb = nbformat.read(src, as_version=4)
client = NotebookClient(nb, timeout=timeout, kernel_name="python3",
                        allow_errors=True, resources={"metadata": {"path": "."}})
t0 = time.time()
client.execute()
nbformat.write(nb, dst)
print(f"\n执行完毕：{src} -> {dst}，用时 {time.time()-t0:.0f} s")

nerr = 0
for i, c in enumerate(nb.cells):
    if c.cell_type != "code":
        continue
    for o in c.get("outputs", []):
        if o.get("output_type") == "error":
            nerr += 1
            print(f"\n===== cell {i} 报错：{o.get('ename')}: {o.get('evalue')}")
            print("\n".join(o.get("traceback", []))[-3000:])
        elif o.get("output_type") == "stream" and o.get("name") == "stderr":
            txt = "".join(o.get("text", ""))
            if "Warning" not in txt and txt.strip():
                print(f"[cell {i} stderr] {txt[:500]}")
print(f"\n报错 cell 数：{nerr}")
