"""Демонстрационный интерфейс расчёта дескрипторов по SMILES."""
import streamlit as st
from tools.rdkit_tools import calculate_descriptors, draw_molecule

st.set_page_config(page_title="Анализ молекулы", layout="wide")
st.title("Анализ молекулы")
st.caption("Первый научный инструмент: расчёт дескрипторов RDKit по введённому SMILES.")
st.info("LLM, агенты, научные базы и остальные сценарии пока не подключены.")
smiles = st.text_input("SMILES", placeholder="Например: CCO")
if st.button("Проанализировать"):
    result = calculate_descriptors(smiles)
    if result.status != "success":
        for error in result.errors:
            st.error(error)
    else:
        st.image(draw_molecule(smiles), caption="Введённая структура")
        st.write("Канонический SMILES:", result.evidence[0].context["canonical_smiles"])
        st.caption("Результаты — расчётные характеристики, а не данные об активности.")
        st.table([{"Показатель": item.method, "Значение": item.value,
                   "Единица": item.unit or "—"} for item in result.evidence])
        st.write("Инструмент: RDKit", result.evidence[0].context["rdkit_version"])
        for limitation in result.limitations:
            st.warning(limitation)
# TODO (IT 1 + IT 2): после выбора LLM подключить workflow и агентные сценарии.
