import streamlit as st

DESCRIPTIONS = {
    "yellow-fever": """
Yellow fever virus (YFV) is a mosquito-borne flavivirus that causes yellow fever, an acute viral hemorrhagic disease found in tropical and subtropical areas of Africa and South America.
It has a single-stranded, positive-sense RNA genome. Outbreaks are historically important and vaccination is an effective preventive measure.
    """,
    "zika": """
Zika virus (ZIKV) is another flavivirus transmitted by *Aedes* mosquitoes.
It caused large outbreaks in the Americas during 2015–2016 and is linked to microcephaly and other neurological disorders.
The virus has a positive-sense RNA genome and is closely related to dengue and yellow fever viruses.
    """,
    "monkeypox": """
Monkeypox virus (MPXV) is a large double-stranded DNA virus belonging to the Orthopoxvirus genus.
It causes monkeypox disease in humans, which is clinically similar to smallpox but generally less severe.
Since 2022, outbreaks have spread internationally, making it an emerging public health concern.
    """,
    "influenza": """
Influenza viruses are segmented negative-sense RNA viruses from the Orthomyxoviridae family.
They cause seasonal flu epidemics and occasional pandemics.
The genome is divided into multiple RNA segments, which enables reassortment and rapid evolution.
    """,
    "sars-cov-2": """
SARS-CoV-2 is the causative agent of the COVID-19 pandemic.
It is a positive-sense single-stranded RNA coronavirus with a ~30 kb genome.
It has generated multiple variants of concern since 2019, leading to extensive global sequencing and surveillance efforts.
    """,
    "rsv-a": """
Respiratory Syncytial Virus (RSV) is a major cause of lower respiratory tract infections in infants and elderly adults.
RSV-A is one of two major antigenic subgroups, with a single-stranded negative-sense RNA genome.
    """,
    "rsv-b": """
Respiratory Syncytial Virus subgroup B (RSV-B) is genetically distinct from RSV-A but causes similar respiratory illness.
It co-circulates with RSV-A and contributes to seasonal epidemics worldwide.
    """,
}

SOURCE_EXPLANATIONS = {
    "ncbi": """
Data is obtained from **NCBI Virus**, using either the CLI or FTP interface.
Sequences and metadata are downloaded directly from the NCBI datasets service.
    """,
    "nextstrain": """
Data is obtained from **Nextstrain open datasets**.
Nextstrain continuously updates curated genomic data for pathogens such as SARS-CoV-2 and RSV.
    """,
    "ftp": """
Data is retrieved from **NCBI FTP servers**, specifically bulk datasets of viral genomes and metadata.
    """
}

REFERENCES = {
    "yellow-fever": [
        {"text": "Yellow Fever - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/yellow-fever"},
        {"text": "Yellow Fever Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_002031"}
    ],
    "zika": [
        {"text": "Zika Virus - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/zika-virus"},
        {"text": "Zika Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_012532"}
    ],
    "monkeypox": [
        {"text": "Monkeypox - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/monkeypox"},
        {"text": "Monkeypox Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_063383"}
    ],
    "influenza": [
        {"text": "Influenza - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/influenza-(seasonal)"},
        {"text": "Influenza Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_026431"}
    ],
    "sars-cov-2": [
        {"text": "COVID-19 - WHO", "link": "https://www.who.int/emergencies/diseases/novel-coronavirus-2019"},
        {"text": "SARS-CoV-2 - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_045512"}
    ],
    "rsv-a": [
        {"text": "RSV - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/rsv"},
        {"text": "RSV-A - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_001781"}
    ],
    "rsv-b": [
        {"text": "RSV - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/rsv"},
        {"text": "RSV-B - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_001782"}
    ]
}

def about(virus):
    st.info(DESCRIPTIONS.get(virus, "Description not available."))

def reference(virus, config):
    st.subheader("Reference Genome")
    ref = config["viruses"][virus].get("reference", {})
    accession = ref.get("accession_id", "N/A")
    length = ref.get("length", "N/A")
    if accession != "N/A":
        st.markdown(
            f"- **Accession ID**: [{accession}](https://www.ncbi.nlm.nih.gov/nuccore/{accession})"
        )
    else:
        st.markdown(f"- **Accession ID**: Not available")
    st.markdown(f"- **Genome Length**: {length:,} bp")

def source(virus, config):
    st.subheader("Data Source")
    src = config["viruses"][virus].get("source", "").lower()
    explanation = SOURCE_EXPLANATIONS.get(src, "Source not documented.")
    st.write(explanation)

def quality_filters(virus, config):
    st.subheader("Preprocessing & Quality Filters")
    filters = config["viruses"][virus]["parameters"].get("metadata_processing", {}).get("filters", [])
    if not filters:
        st.write("No explicit filters defined for this virus.")
        return
    st.write("The following filters are applied to raw metadata before analysis:")
    for rule in filters:
        column = rule.get("column")
        op = rule.get("operator")
        val = rule.get("value", "")
        st.markdown(f"- **{column}** {op} {val}".strip())

def haplocov_parameters(virus, config):
    st.subheader("HaploCov Parameters")
    haplo = config["viruses"][virus]["parameters"].get("haplocov", {})
    if not haplo:
        st.write("No HaploCov parameters defined for this virus.")
        return
    st.write("HaploCov is a tool that designates viral haplotypes based on clustering of genomes. Parameters used:")
    st.markdown(f"- **Distance threshold (dist):** {haplo.get('dist')}")
    st.markdown(f"- **Minimum cluster size (size):** {haplo.get('size')}")
    st.markdown(f"- **Designation mode:** {haplo.get('designation_mode')}")

def dataset(virus, df):
    st.subheader("Dataset Overview")
    if df is None or df.empty:
        st.warning("No dataset available.")
        return
    st.write("The dataset contains the following columns:")
    st.write(list(df.columns))
    st.write(f"Total records: {len(df):,}")

def references(virus):
    st.subheader("References & Resources")
    refs = REFERENCES.get(virus, [])
    if not refs:
        st.write("No references available.")
        return
    for ref in refs:
        st.markdown(f"- [{ref['text']}]({ref['link']})")

def describe(virus, config, df):
    about(virus)
    st.markdown("---")

    reference(virus, config)
    st.markdown("---")

    source(virus, config)
    st.markdown("---")

    quality_filters(virus, config)
    st.markdown("---")

    haplocov_parameters(virus, config)
    st.markdown("---")

    dataset(virus, df)
    st.markdown("---")

    references(virus)
