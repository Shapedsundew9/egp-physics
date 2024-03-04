"""FGC Group class for EGP Physics."""
from __future__ import annotations
from typing import cast
from random import randint
from egp_types.aGC import aGC
from egp_types.mermaid_charts import MERMAID_IGRAPH_CLASS_DEF_STR
from .egp_typing import InsertRow


class insertion_work:
    """FGC Group class for EGP Physics.

    Stores the results of recursive insertion in a structure that is easy to resolve
    the references in.
    """

    def __init__(
        self,
        tgc: aGC,
        igc: aGC,
        rgc: insertion_work | aGC | None,
        fgc: insertion_work | aGC | None,
        gca: bool | None = None,
        above_row: InsertRow = "A",
    ) -> None:
        self.tgc: aGC = tgc
        self.igc: aGC = igc
        self.rgc: insertion_work | aGC | None = rgc
        self.fgc: insertion_work | aGC | None = fgc
        self.gca: bool | None = gca
        self.above_row: InsertRow = above_row
        self._rgc: insertion_work | None = None
        self._fgc: insertion_work | None = None

    def __repr__(self) -> str:
        """Return a non-recursive representation of the insertion_work."""
        ret_str = f"insertion_work.tgc\n{self.tgc!r}"
        ret_str += f"\ninsertion_work.igc\n{self.igc!r}"
        if isinstance(self.rgc, insertion_work):
            ret_str += "\ninsertion_work.rgc\nUNRESOLVED\n"
        elif self.rgc is not None:
            ret_str += f"\ninsertion_work.rgc\n{self.rgc!r}"
        if isinstance(self.fgc, insertion_work):
            ret_str += "\ninsertion_work.fgc\nUNRESOLVED\n"
        elif self.fgc is not None:
            ret_str += f"\ninsertion_work.fgc\n{self.fgc!r}"
        if self.gca is not None:
            ret_str += f"\ninsertion_work.gca\n{self.gca!r}"
        ret_str += f"\ninsertion_work.above_row\n{self.above_row!r}"
        return ret_str

    def mermaid_chart(self) -> str:
        """Return a mermaid chart of the insertion_work.
        Work is recursively embedded in subgraphs to give a holistic view.
        Chart can be viewed at https://mermaid.live/
        """
        subgraph_strs: list[str]
        link_strs: list[str]
        subgraph_strs, link_strs = self.mermaid_chart_body()
        ret_str: str = "\n".join(subgraph_strs)
        ret_str += "\n\n" + "\n".join(link_strs)
        ret_str += "\n\n" + "\n".join(self.tgc["gc_graph"].i_graph.mermaid_style_str(link_strs))
        uids: set[int] = {int(link_str[-4:], 16) for link_str in link_strs}
        ret_str += "\n\n" + MERMAID_IGRAPH_CLASS_DEF_STR
        ret_str += "\n\n" + "\n".join(("\n".join(self.igc["gc_graph"].i_graph.mermaid_class_str(uid)) for uid in uids))
        return ret_str

    def mermaid_chart_body(self, top: bool = True) -> tuple[list[str], list[str]]:
        """Return a mermaid chart of the insertion_work.
        If RGC or FGC are themselves defined by insertion work then they are recursively embedded in subgraphs to give a holistic view.
        """
        tgc_strs: tuple[list[str], list[str]] = self.tgc["gc_graph"].i_graph.mermaid_embedded_str(randint(0, 2**16 - 1))
        igc_strs: tuple[list[str], list[str]] = self.igc["gc_graph"].i_graph.mermaid_embedded_str(randint(0, 2**16 - 1))
        rgc_strs: tuple[list[str], list[str]] = ([], [])
        if isinstance(self._rgc, insertion_work):
            rgc_strs = self._rgc.mermaid_chart_body(False)
        elif self.rgc is not None:
            rgc_strs = cast(aGC, self.rgc)["gc_graph"].i_graph.mermaid_embedded_str(randint(0, 2**16 - 1))
        fgc_strs: tuple[list[str], list[str]] = ([], [])
        if isinstance(self._fgc, insertion_work):
            fgc_strs = self._fgc.mermaid_chart_body(False)
        elif self.fgc is not None:
            fgc_strs = cast(aGC, self.fgc)["gc_graph"].i_graph.mermaid_embedded_str(randint(0, 2**16 - 1))

        # Contain each GC in a uniquely identifable subgraph
        uid: int = randint(0, 2**16 - 1)
        ret_strs: tuple[list[str], list[str]] = ([], [])
        if top:
            ret_strs[0].append("flowchart TB")
        ret_strs[0].append(f'\tsubgraph W{uid:04x}["Work"]')
        for name, strs in [("TGC", tgc_strs), ("IGC", igc_strs), ("RGC", rgc_strs), ("FGC", fgc_strs)]:
            for i in range(len(strs[0])):
                strs[0][i] = f"\t\t{strs[0][i]}"
            if strs[0]:
                strs[0].insert(0, f'\t\tsubgraph {name}{uid:04x}["{name}"]')
                strs[0].append("\t\tend")
            ret_strs[0].extend(strs[0])
            ret_strs[1].extend(strs[1])
        ret_strs[0].append("\tend")
        return ret_strs

    def resolve(self) -> dict[int, aGC]:
        """Resolve the insertion work into a dictionary of GCs.

        Returned dictionary has the structure
        reference: GC
        with the order [tgc, igc, _rgc_, _fgc_] where
        _rgc_ is either an rgc or [tgc, igc, _rgc_, _fgc_]
        _fgc_ may be None, an fgc or [tgc, igc, _rgc_, _fgc_]
        """
        retval: dict[int, aGC] = {self.tgc["ref"]: self.tgc, self.igc["ref"]: self.igc}
        if isinstance(self.rgc, insertion_work):
            retval.update(self.rgc.resolve())
            self._rgc = self.rgc
            self.rgc = cast(aGC, self.rgc.rgc)
        if self.rgc is not None:
            retval[self.rgc["ref"]] = self.rgc
            if isinstance(self.fgc, insertion_work):
                retval.update(self.fgc.resolve())
                assert self.gca is not None
                self._fgc = self.fgc
                self.fgc = cast(aGC, self.fgc.rgc)
                self.rgc[("gcb_ref", "gca_ref")[self.gca]] = self.fgc["ref"]
            elif self.fgc is not None:
                retval[self.fgc["ref"]] = self.fgc
        return retval
