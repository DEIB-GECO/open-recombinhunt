from pathlib import Path
import sys

import pandas as pd

SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))
from src.utils.constants import *

import streamlit as st
import plotly.express as px

INFO = {
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

DESCRIPTIONS = {
    "yellow-fever": """
Yellow fever is an epidemic-prone mosquito-borne vaccine preventable disease that is transmitted to humans by the bites of infected mosquitoes. Yellow fever is caused by an arbovirus (a virus transmitted by vectors such mosquitoes, ticks or other arthropods) transmitted to humans by the bites of infected Aedes and Haemagogus mosquitoes.

These day-biting mosquitoes breed around houses (domestic), in forests or jungles (sylvatic), or in both habitats (semi-domestic). Yellow fever is a high-impact high-threat disease, with risk of international spread, which represents a potential threat to global health security.

The incubation period for yellow fever is 3 to 6 days. Many people do not experience symptoms. Common symptoms include fever, muscle pain, headache, loss of appetite, nausea or vomiting. In most cases, symptoms disappear after 3 to 4 days.

A small percentage of patients enter a second, more toxic phase within 24 hours of recovering from initial symptoms. High fever returns and several body systems are affected, usually the liver and the kidneys. In this phase, people are likely to develop jaundice (yellowing of the skin and eyes, hence the name yellow fever), dark urine, and abdominal pain with vomiting. Bleeding can occur from the mouth, nose, eyes, or stomach. Half of the patients who enter the toxic phase die within 7–10 days.

There is no specific anti-viral drug for yellow fever. Patients should rest, stay hydrated and seek medical advice. Depending on the clinical manifestations and other circumstances, patients may be sent home, be referred for in-hospital management, or require emergency treatment and urgent referral. Treatment for dehydration, liver and kidney failure, and fever improves outcomes. Associated bacterial infections can be treated with antibiotics.
    """,
    "zika": """
Zika virus is most commonly spread to people by the bite of an infected Aedes species mosquito. It can also be spread through sex from a person who is infected with Zika virus to their sexual partner(s).

Zika virus can be passed from a pregnant woman to her fetus. Infection during pregnancy can cause certain birth defects.

Many people infected with Zika will not have symptoms or will only have mild symptoms. The most common symptoms are fever, rash, headache, joint and muscle pain, and red eyes.

Zika virus typically occurs in tropical and subtropical areas of Africa, the Americas, Southern Asia, and Western Pacific.

There is currently no vaccine to prevent or medicine to treat Zika.

Travelers can protect themselves by preventing mosquito bites. When traveling to countries with Zika virus:

- Use EPA-registered insect repellent.
- Wear long-sleeved shirts and pants.
- Stay in places with air conditioning or that use window and door screens.
    """,
    "monkeypox": """
Mpox is an infectious disease that can cause a painful rash, enlarged lymph nodes, fever, headache, muscle ache, back pain and low energy. Most people fully recover, but some get very sick. 

Mpox is caused by the monkeypox virus (MPXV). It is an enveloped double-stranded DNA virus of the Orthopoxvirus genus in the Poxviridae family, which includes variola, cowpox, vaccinia and other viruses. There are two distinct clades of the virus: clade I (with subclades Ia and Ib) and clade II (with subclades IIa and IIb).
A global outbreak of clade IIb began in 2022 and continues to this day, including in some African countries. There are also growing outbreaks of clades Ia and Ib affecting the Democratic Republic of the Congo and other countries in Africa. As of August 2024, clade Ib has also been detected beyond Africa.
The natural reservoir of the virus is unknown, but various small mammals such as squirrels and monkeys are susceptible. 

Mpox spreads from person to person mainly through close contact with someone who has mpox, including members of a household. Close contact includes skin-to-skin (such as touching or sex) and mouth-to-mouth or mouth-to-skin contact (such as kissing), and it can also include being face-to-face with someone who has mpox (such as talking or breathing close to one another, which can generate infectious respiratory particles).
People with multiple sexual partners are at higher risk of acquiring mpox. 
People can also contract mpox from contaminated objects such as clothing or linen, through needle injuries in health care, or in community settings such as tattoo parlours. 

During pregnancy or birth, the virus may be passed to the baby. Contracting mpox during pregnancy can be dangerous for the fetus or newborn infant and can lead to loss of the pregnancy, stillbirth, death of the newborn, or complications for the parent.
Animal-to-human transmission of mpox occurs from infected animals to humans from bites or scratches, or during activities such as hunting, skinning, trapping, cooking, playing with carcasses or eating animals. The animal reservoir of the monkeypox virus remains unknown and further studies are underway. 
    """,
    "influenza": """
H5N1 is one of several influenza viruses that causes a highly infectious respiratory disease in birds called avian influenza (or "bird flu"). Infections in mammals, including humans, have also been documented.

H5N1 influenza virus infection can cause a range of diseases in humans, from mild to severe and in some cases, it can even be fatal. Symptoms reported have primarily been respiratory, but conjunctivitis and other non-respiratory symptoms have also been reported. There have also been a few detections of A(H5N1) virus in persons who were exposed to infected animals or their environments but who did not show any symptoms.  

The goose/Guangdong-lineage of H5N1 avian influenza viruses first emerged in 1996 and has been causing outbreaks in birds since then. Since 2020, a variant of these viruses has led to an unprecedented number of deaths in wild birds and poultry in many countries. First affecting Africa, Asia and Europe, in 2021, the virus spread to North America, and in 2022, to Central and South America. From 2021 to 2022, Europe and North America observed their largest and most extended epidemic of avian influenza with unusual persistence of the virus in wild bird populations.

Since 2022, there have been increasing reports of deadly outbreaks among mammals also caused by influenza A(H5) – including influenza A(H5N1) – viruses. There are likely to be more outbreaks that have not been detected or reported. Both land and sea mammals have been affected, including outbreaks in farmed fur animals, seals, sea lions, and detections in other wild and domestic animals such as foxes, bears, otters, raccoons, cats, dogs, cows, goats and others.
    """,
    "sars-cov-2": """
COVID-19 is the disease caused by the SARS-CoV-2 coronavirus. It usually spreads between people in close contact. COVID-19 vaccines provide strong protection against severe illness and death. Although a person can still get COVID-19 after vaccination, they are more likely to have mild or no symptoms. Anyone can get sick with COVID-19 and become seriously ill or die, but most people will recover without treatment. People over age 60 and those with existing medical conditions have a higher risk of getting seriously ill. These conditions include high blood pressure, diabetes, obesity, immunosuppression including HIV, cancer and pregnancy. Unvaccinated people also have a higher risk of severe symptoms. 

People may experience different symptoms from COVID-19. Symptoms usually begin 5–6 days after exposure and last 1–14 days.

The most common symptoms are:

- fever
- chills
- sore throat

Less common symptoms are:

- muscle aches and heavy arms or legs
- severe fatigue or tiredness
- runny or blocked nose, or sneezing
- headache
- sore eyes
- dizziness
- new and persistent cough
- tight chest or chest pain
- shortness of breath
- hoarse voice
- numbness or tingling
- appetite loss, nausea, vomiting, abdominal pain or diarrhoea
- loss or change of sense of taste or smell
- difficulty sleeping.

People with the following symptoms should seek immediate medical attention:

- difficulty breathing, especially at rest, or unable to speak in sentences
- confusion
- drowsiness or loss of consciousness
- persistent pain or pressure in the chest
- skin being cold or clammy, or turning pale or a bluish colour
- loss of speech or movement.

People who have pre-existing health problems are at higher risk when they have COVID-19; they should seek medical help early if worried about their condition. These include people taking immunosuppressive medication; those with chronic heart, lung, liver or rheumatological problems; those with HIV, diabetes, cancer. obesity or dementia.
People with severe disease and those needing hospital treatment should receive treatment as soon as possible. The consequences of severe COVID-19 include death, respiratory failure, sepsis, thromboembolism (blood clots), and multiorgan failure, including injury of the heart, liver or kidneys.
In rare situations, children can develop a severe inflammatory syndrome a few weeks after infection. 
Some people who have had COVID-19, whether they have needed hospitalization or not, continue to experience symptoms. These long-term effects are called long COVID (or post COVID-19 condition). The most common symptoms associated with long COVID include fatigue, breathlessness and cognitive dysfunction (for example, confusion, forgetfulness, or a lack of mental focus or clarity). Long COVID can affect a person’s ability to perform daily activities such as work or household chores.  

    """,
    "rsv-a": """
Respiratory syncytial virus (RSV) belongs to the genus Orthopneumovirus within the family Pneumoviridae and order Mononegavirales. Members of this genus include human RSV, bovine RSV and murine pneumonia virus. There are two major antigenic subtypes of human RSV (A and B) determined largely by antigenic drift and duplications in RSV-G sequences, but accompanied by genome-wide sequence divergence, including within RSV-F.

Human RSV is a globally prevalent cause of lower respiratory tract infection in all age groups. In infants and young children, the first infection may cause severe bronchiolitis that can sometimes be fatal. In older children and adults without comorbidities, repeated upper respiratory tract infections are common and range from subclinical infection to symptomatic upper respiratory tract disease. In addition to the pediatric burden of disease, RSV is increasingly being recognized as an important pathogen in older adults, with infection leading to an increase in hospitalization rates among those aged 65 years and over, and to increased mortality rates among the frail elderly that approach the rates seen with influenza. The risk of severe disease in adults is increased by the presence of underlying chronic pulmonary disease, circulatory conditions and functional disability, and is associated with higher viral loads. RSV is also a nosocomial threat both to young infants and among immunocompromised and vulnerable individuals. High mortality rates have been observed in those infected with RSV following bone marrow or lung transplantation.
    """,
    "rsv-b": """
Respiratory syncytial virus (RSV) belongs to the genus Orthopneumovirus within the family Pneumoviridae and order Mononegavirales. Members of this genus include human RSV, bovine RSV and murine pneumonia virus. There are two major antigenic subtypes of human RSV (A and B) determined largely by antigenic drift and duplications in RSV-G sequences, but accompanied by genome-wide sequence divergence, including within RSV-F.

Human RSV is a globally prevalent cause of lower respiratory tract infection in all age groups. In infants and young children, the first infection may cause severe bronchiolitis that can sometimes be fatal. In older children and adults without comorbidities, repeated upper respiratory tract infections are common and range from subclinical infection to symptomatic upper respiratory tract disease. In addition to the pediatric burden of disease, RSV is increasingly being recognized as an important pathogen in older adults, with infection leading to an increase in hospitalization rates among those aged 65 years and over, and to increased mortality rates among the frail elderly that approach the rates seen with influenza. The risk of severe disease in adults is increased by the presence of underlying chronic pulmonary disease, circulatory conditions and functional disability, and is associated with higher viral loads. RSV is also a nosocomial threat both to young infants and among immunocompromised and vulnerable individuals. High mortality rates have been observed in those infected with RSV following bone marrow or lung transplantation.
    """,
}

REFERENCES = {
    "yellow-fever": [
        {"text": "Yellow Fever - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/yellow-fever"},
        {"text": "Yellow Fever Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_002031"}
    ],
    "zika": [
        {"text": "Zika Virus - CDC", "link": "https://www.cdc.gov/zika/about/index.html"},
        {"text": "Zika Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_012532"}
    ],
    "monkeypox": [
        {"text": "Mpox - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/mpox"},
        {"text": "Monkeypox Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_063383"}
    ],
    "influenza": [
        {"text": "Influenza H5N1 - WHO", "link": "https://www.who.int/news-room/questions-and-answers/item/influenza-h5n1"},
        {"text": "Influenza Virus - NCBI", "link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_026431"}
    ],
    "sars-cov-2": [
        {"text": "COVID-19 - WHO", "link": "https://www.who.int/news-room/fact-sheets/detail/coronavirus-disease-(covid-19)"},
    ],
    "rsv-a": [
        {"text": "RSV - WHO", "link": "https://www.who.int/teams/health-product-policy-and-standards/standards-and-specifications/norms-and-standards/vaccine-standardization/respiratory-syncytial-virus-disease"},
    ],
    "rsv-b": [
        {"text": "RSV - WHO", "link": "https://www.who.int/teams/health-product-policy-and-standards/standards-and-specifications/norms-and-standards/vaccine-standardization/respiratory-syncytial-virus-disease"},
    ],
    "haplocov": [
        {"text": "HaploCoV - GitHub", "link": "https://github.com/matteo14c/HaploCoV"},
        {"text": "HaploCoV - Publication", "link": "https://www.nature.com/articles/s42003-023-04784-4"}
    ]
}

def about(virus):
    st.info(INFO.get(virus))
    st.write(DESCRIPTIONS.get(virus))

def reference(virus, config):
    st.subheader("Reference Genome")
    st.write("The reference genome serves as a standard for aligning and comparing viral sequences. It is typically a well-characterized isolate that represents the species or strain of interest.")
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

    # add button to download reference genome
    # reference.fasta file is generated in data/processed/{virus}/reference.fasta
    # get path from config: paths.processed_data
    processed_path = config[PATHS].get(PROCESSED_DATA)
    ref_path = f"{processed_path}/{virus}/reference.fasta"

    # if file does not exist, skip
    try:
        with open(ref_path, "r") as f:
            fasta_data = f.read()
        st.download_button(
            label="Download Reference Genome (FASTA)",
            data=fasta_data,
            file_name=f"{virus}_reference_{accession}.fasta",
            mime="text/plain"
        )
    except FileNotFoundError:
        pass

SOURCE_LINKS = {
    "yellow-fever": "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Yellow%20fever%20virus,%20taxid:11089",
    "zika": "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Zika%20virus,%20taxid:64320",
    "monkeypox": "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Monkeypox%20virus,%20taxid:10244",
    "rsv-a": "https://nextstrain.org/rsv/a/genome/all-time",
    "rsv-b": "https://nextstrain.org/rsv/b/genome/all-time",
    "sars-cov-2": "https://nextstrain.org/ncov/open/global/6m",
}

def source(virus, config):
    st.subheader("Data Source")
    src = config[VIRUSES][virus].get(SOURCE).lower()

    if src == NCBI:
        # write about NCBI, provide link to NCBI Virus
        st.markdown("Data is sourced from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), a comprehensive database of viral sequences and related information.")

        # mention the taxon id used to fetch data and provide link to NCBI Taxonomy
        taxon_id = config[VIRUSES][virus].get(TAXON_ID)
        if taxon_id:
            st.markdown(f"- NCBI Taxonomy Browser for Taxon ID used for data retrieval: [{taxon_id}]({SOURCE_LINKS.get(virus)})")

        # mention the CLI tool used to fetch data
        cli_sequences = config[NCBI_CLI].get(SEQUENCES)
        cli_reference = config[NCBI_CLI].get(REFERENCE)

        # write it like a code block
        st.markdown("The following NCBI CLI commands are used to fetch data:")
        st.code(f"""
        # Command to download the main dataset.
        {cli_sequences[0].format(taxon_id=taxon_id, virus_name=virus)}
        # Command to generate metadata from the downloaded report.
        {cli_sequences[1].format(taxon_id=taxon_id, virus_name=virus)}

        # Command to download the reference genome.
        {cli_reference[0].format(accession_id=config[VIRUSES][virus].get(REFERENCE).get(ACCESSION_ID), virus_name=virus)}
        """, language="bash")

    elif src == NEXTSTRAIN:
        # write about Nextstrain, provide link to Nextstrain
        st.markdown("Data is sourced from [Nextstrain](https://nextstrain.org/), an open-source project that provides real-time tracking of pathogen evolution.")

        # provide links to the specific virus page on Nextstrain
        nextstrain_url = SOURCE_LINKS.get(virus)
        if nextstrain_url:
            st.markdown(f"- Nextstrain page for {virus}: [{nextstrain_url}]({nextstrain_url})")

        # mention the URL to fetch data
        url_sequences = config[NEXTSTRAIN_URL].get(virus).get(SEQUENCES, None)
        url_metadata = config[NEXTSTRAIN_URL].get(virus).get(METADATA, None)

        if url_metadata:
            st.markdown(f"- URL to download metadata: {url_metadata}")
        if not url_sequences:
            st.warning("No sequences are downloaded for sars-cov-2 as its metadata already contains the list of mutations per sequence. Hence, only metadata is downloaded since we will not perform sequence-level analysis by HaploCoV.")
        if url_sequences:
            st.markdown(f"- URL to download sequences: {url_sequences}")

            # even if the metadata and sequences are from URLs, we still use the ncbi.cli tool to download reference sequence
            st.info("Note: The reference genome is still fetched using the NCBI CLI tool as described below.")
            cli_reference = config[NCBI_CLI].get(REFERENCE)
            st.code(f"""
            # Command to download the reference genome.
            {cli_reference[0].format(accession_id=config[VIRUSES][virus].get(REFERENCE).get(ACCESSION_ID), virus_name=virus)}
            """, language="bash")

    elif src == FTP:
        # write about FTP, provide link to NCBI FTP
        st.markdown("Data is sourced from [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/), a repository for various biological data including genomic sequences.")

        # mention that we used it because it allowed us to fetch the data that we require.
        # we needed Influenza H5N1 sequences of HA segment, that are collected from North America only. NCBI Virus did not allow us to filter accordingly.
        st.info("Note: The FTP source was chosen because it allowed for specific filtering of sequences, such as obtaining Influenza H5N1 HA segment sequences collected from North America, which was not feasible through NCBI Virus.")

        # in this case we used a combination of different sources to fetch data,
        # actually only metadata is fetched from FTP, sequences are fetched from NCBI Virus using the accession ids after filtering the metadata.
        # reference genome is fetched using NCBI CLI tool as well.
        st.markdown("In this case, metadata is fetched from the FTP source, while sequences are obtained from NCBI Virus using the filtered accession IDs. The reference genome is also fetched using the NCBI CLI tool as described below.")
        ftp_url = config[FTP_URL].get(virus).get(METADATA, None)
        if ftp_url:
            st.markdown(f"- URL to download metadata: {ftp_url}")

        cli_sequences = config[FTP_CLI].get(SEQUENCES)
        cli_reference = config[NCBI_CLI].get(REFERENCE)
        st.markdown("The following commands are used to fetch remaining data:")
        st.code(f"""
        # Command to download the sequences from NCBI.
        {cli_sequences[0].format(virus_name=virus)}
        
        # Command to download the reference genome.
        {cli_reference[0].format(accession_id=config[VIRUSES][virus].get(REFERENCE).get(ACCESSION_ID), virus_name=virus)}
        """, language="bash") 

def quality_filters(virus, config):
    
    # after data is fetched, we apply quality filters to remove low-quality sequences.
    # here you can find the list of filters applied to {this virus}.
    # filters are stored in config.yaml for easy access.
    st.subheader("Quality Filters")

    st.write("After data is fetched, quality filters are applied to remove low-quality sequences.")

    if virus == "yellow-fever":
        # Isolate Collection date is not NA
        # Accession is not NA
        # Length is not NA and >= 9775. Note that the reference genome length is 10862, but we set the threshold to 90% of it.
        st.markdown("The following quality filters are applied to the yellow-fever dataset:")
        st.markdown("- Isolate Collection date is not NA")
        st.markdown("- Accession is not NA")
        st.markdown("- Length is not NA and >= 9,775 bp (90% of reference genome length)")

    elif virus == "zika":
        """
        filters:
          # Rule 1: Discard rows where 'Isolate Collection date' is missing.
          - column: "Isolate Collection date"
            operator: "notna"

          # Rule 2: Discard rows where 'Accession' is missing.
          - column: "Accession"
            operator: "notna"

          # Rule 3: Discard rows where 'Length' is missing.
          - column: "Length"
            operator: "notna"

          # Rule 4: Discard rows where sequence 'Length' is less than 90% of the reference length.
          - column: "Length"
            operator: ">="
            value: 9727"""
        st.markdown("The following quality filters are applied to the zika dataset:")
        st.markdown("- Isolate Collection date is not NA")
        st.markdown("- Accession is not NA")
        st.markdown("- Length is not NA and >= 9,727 bp (90% of reference genome length)")

    elif virus == "monkeypox":
        """
        filters:
          # Rule 1: Discard rows where 'Isolate Collection date' is missing.
          - column: "Isolate Collection date"
            operator: "notna"

          # Rule 2: Discard rows where 'Accession' is missing.
          - column: "Accession"
            operator: "notna"

          # Rule 3: Discard rows where 'Length' is missing.
          - column: "Length"
            operator: "notna"

          # Rule 4: Discard rows where sequence 'Length' is less than 180000.
          - column: "Length"
            operator: ">="
            value: 177488"""
        st.markdown("The following quality filters are applied to the monkeypox dataset:")
        st.markdown("- Isolate Collection date is not NA")
        st.markdown("- Accession is not NA")
        st.markdown("- Length is not NA and >= 177,488 bp (90% of reference genome length)")

    elif virus == "influenza":
        """
            filtered_df = df[
            (df["Species"] == "Alphainfluenzavirus influenzae") &
            (df["Genotype"].str.contains("H5N1|h5n1", na=False)) &
            (df["Segment"].str.contains("HA|segment 4|Segment 4|4", na=False)) &
            (
                (df["Geo_Location"].str.contains("USA|Canada|Mexico", na=False) & ~df["Geo_Location"].str.contains(":", na=False)) |
                (df["Geo_Location"].str.contains("USA:|Canada:|Mexico:", na=False) & df["Geo_Location"].str.contains(":", na=False))
            ) &
            (df["Length"] > 1672) &
            (df["Collection_Date"].notna()) &
            (df["Release_Date"].notna())
        ].copy()"""
        st.markdown("The following quality filters are applied to the influenza dataset:")
        st.markdown("- Species is Alphainfluenzavirus influenzae")
        st.markdown("- Genotype contains H5N1")
        st.markdown("- Segment contains HA or segment 4")
        st.markdown("- Geo_Location is from North America (USA, Canada, Mexico)")
        st.markdown("- Length is not NA and > 1,672 bp (95% of reference genome length)")
    
    elif virus == "sars-cov-2":
        """
        filters:
          - {column: "missing_data", operator: "<", value: 589} # 2% of 29903
          - {column: "coverage", operator: ">=", value: 0.99}
          - {column: "virus", operator: "==", value: "ncov"}
          - {column: "virus", operator: "notna"}
          - {column: "length", operator: "notna"}
          - {column: "date_submitted", operator: "notna"}
          - {column: "QC_overall_status", operator: "notna"}
          - {column: "QC_overall_status", operator: "!=", value: "bad"}
          - {column: "QC_missing_data", operator: "==", value: "good"}
          - {column: "QC_frame_shifts", operator: "==", value: "good"}
          - {column: "QC_stop_codons", operator: "==", value: "good"}
          - {column: "QC_mixed_sites", operator: "==", value: "good"}"""
        st.markdown("The following quality filters are applied to the sars-cov-2 dataset:")
        st.markdown("- Missing data < 589 bases (2% of 29,903 bp -reference genome length-)")
        st.markdown("- Coverage >= 99%")
        st.markdown("- Virus is ncov")
        st.markdown("- Virus is not NA")
        st.markdown("- Length is not NA")
        st.markdown("- Date submitted is not NA")
        st.markdown("- QC overall status is not NA and not 'bad'")
        st.markdown("- QC missing data is 'good'")
        st.markdown("- QC frame shifts is 'good'")
        st.markdown("- QC stop codons is 'good'")
        st.markdown("- QC mixed sites is 'good'")

    elif virus in ["rsv-a", "rsv-b"]:
        """
        filters:
          - { column: "accession", operator: "notna" }
          - { column: "date", operator: "notna" }
          - { column: "qc.overallStatus", operator: "==", value: "good" }
          - { column: "missing_data", operator: "<=", value: 3 }"""
        st.markdown(f"The following quality filters are applied to the {virus} dataset:")
        st.markdown("- Accession is not NA")
        st.markdown("- Date is not NA")
        st.markdown("- QC overall status is 'good'")
        st.markdown("- Missing data <= 3 bases -decided by looking at the distribution of missing data in the dataset-")
def haplocov():
    st.subheader("What is HaploCoV?")

    # haplocov is basically a tool for inferring viral haplotypes from sequence data.
    # it allows us to assign lineages to sequences based on their genetic similarity.
    # and also identify the mutations that each sequence has faced.

    st.markdown("HaploCov is a novel software framework designed for the unsupervised classification and rapid detection of emerging viral variants. For the OpenRecombinHunt pipeline, HaploCov serves two primary functions for viruses that lack a pre-calculated list of mutations:")
    st.markdown("1. **Mutation Calling**")
    st.markdown("First, HaploCov identifies the specific mutations present in each viral sequence. It takes the processed FASTA sequences and the reference sequence as input and performs a genome alignment for each sequence using the nucmer program. From this alignment, it extracts all genetic variants (substitutions, insertions, and deletions) and produces a comprehensive list of mutations for each genome. This provides the foundational data required for all subsequent analysis.")
    st.markdown("2. **Lineage Designation (Haplogroup Formation)**")
    st.markdown("The second, and more complex, function of HaploCov is to assign each sequence to a meaningful lineage or 'haplogroup' (HG). As described by Chiara et al. (2023), this is achieved through agglomerative hierarchical clustering of phenetic profiles. In simple terms, the tool groups sequences based on their shared patterns of high-frequency mutations. This allows the pipeline to handle two distinct scenarios:")
    st.markdown("- **De Novo Clustering:** For viruses without a pre-existing, reliable classification, HaploCov builds a new classification system from the ground up by clustering the sequences into novel haplogroups.")
    st.markdown("- **Augmentation:** For viruses that already have a baseline classification (like RSV), HaploCov uses the pre-assigned lineages as a starting point and applies the same clustering logic to identify potential sub-clusters, thus refining the existing nomenclature.")
    st.markdown("This process ensures that every virus analyzed in the pipeline has a consistent and meaningful lineage assignment before being passed to the RecombinHunt tool.")
    st.markdown("For more information, please refer to the [HaploCoV GitHub repository](https://github.com/matteo14c/HaploCoV) and the [original publication](https://www.nature.com/articles/s42003-023-04784-4)")

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

def dataset_from_df(virus, df: pd.DataFrame):
    st.subheader("Dataset Overview")
    
    # total number of records in the final dataset
    total_records = len(df)
    st.markdown(f"- **Total Records:** {total_records:,}")

    # number of unique countries in the dataset
    country_col = "Location" if virus == "sars-cov-2" else "country"
    if country_col in df.columns:
        unique_countries = df[country_col].nunique()
        st.markdown(f"- **Unique Countries:** {unique_countries}")

        with st.expander("Country Distribution", expanded=False):
            country_counts = df[country_col].value_counts().reset_index()
            country_counts.columns = ["Country", "Count"]
            fig = px.bar(country_counts, x="Country", y="Count", title="Number of Sequences per Country (Log Scale)" if virus == "sars-cov-2" else "Number of Sequences per Country", log_y=True if virus == "sars-cov-2" else False)
            st.plotly_chart(fig, use_container_width=True)

    # number of unique lineages in the dataset, both named pangoLin
    lineage_col = "pangoLin" 
    if lineage_col in df.columns:
        unique_lineages = df[lineage_col].nunique()
        st.markdown(f"- **Unique Lineages:** {unique_lineages}")

        with st.expander("Lineage Distribution", expanded=False):
            # display it as a table, not a bar plot
            lineage_counts = df[lineage_col].value_counts().reset_index()
            lineage_counts.columns = ["Lineage", "Count"]
            st.dataframe(lineage_counts, hide_index=True, use_container_width=False)

    if virus not in ["sars-cov-2"]:
        # number of lineages assigned by HaploCov
        # count lineages in the lineage_col that contains "NmC"

        # number of lineages before HaploCov is run
        # total number of unique values in lineage_col - number of unique values that contain "NmC"
        if lineage_col in df.columns:
            total_lineages = df[lineage_col].nunique()
            haplocov_lineages = df[df[lineage_col].str.contains("NmC", na=False)][lineage_col]
            existing_lineages = df[~df[lineage_col].str.contains("NmC", na=False)][lineage_col]

            haplocov_lineages_count = haplocov_lineages.nunique()
            existing_lineages_count = existing_lineages.nunique()

            if existing_lineages_count == 1:
                # there is no existing nomenclature for this virus, A.1 assigned to all sequences before HaploCoV is run.
                # so we can say that HaploCov created all the lineages from scratch.
                st.info("HaploCov created all lineages from scratch as there was no existing nomenclature for this virus. \nA.1 was assigned to all sequences before HaploCoV was run.")

                st.markdown(f"- **Lineages assigned by HaploCoV:** {haplocov_lineages_count + existing_lineages_count}")
                with st.expander("HaploCoV-assigned lineages:", expanded=False):
                    all_lineages_df = df[lineage_col].value_counts().reset_index()
                    all_lineages_df.columns = ["Lineage", "Count"]
                    st.dataframe(all_lineages_df, hide_index=True, use_container_width=False)

            else:
                st.markdown(f"- **Lineages before HaploCoV:** {existing_lineages_count}")
                with st.expander("Existing nomenclature:", expanded=False):
                    existing_lineages_df = existing_lineages.value_counts().reset_index()
                    existing_lineages_df.columns = ["Lineage", "Count"]
                    st.dataframe(existing_lineages_df, hide_index=True, use_container_width=False)

                st.markdown(f"- **Lineages assigned by HaploCoV:** {haplocov_lineages_count}")
                with st.expander("HaploCoV-assigned lineages:", expanded=False):
                    haplocov_lineages_df = haplocov_lineages.value_counts().reset_index()
                    haplocov_lineages_df.columns = ["Lineage", "Count"]
                    st.dataframe(haplocov_lineages_df, hide_index=True, use_container_width=False)

def dataset_from_stats(virus, stats: dict):
    st.subheader("Dataset Overview")
    
    # total number of records in the final dataset
    total_records = stats.get("total_records")
    st.markdown(f"- **Total Records:** {total_records:,}")

    min_collection_date = stats.get("min_collection_date", "N/A")
    max_collection_date = stats.get("max_collection_date", "N/A")
    st.markdown(f"- **Collection Date Range:** {min_collection_date} to {max_collection_date}")

    # number of unique countries in the dataset
    unique_countries = stats.get("unique_countries", 0)
    st.markdown(f"- **Unique Countries:** {unique_countries}")

    with st.expander("Country Distribution", expanded=False):
        country_counts = stats.get("country_distribution", {})
        if country_counts is not None and not country_counts.empty:
            country_counts.columns = ["Country", "Count"]
            fig = px.bar(country_counts, x="Country", y="Count", title="Number of Sequences per Country (Log Scale)", log_y=True)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.write("No country data available.")

    # number of unique lineages in the dataset, both named pangoLin
    unique_lineages = stats.get("unique_lineages", 0)
    st.markdown(f"- **Unique Lineages:** {unique_lineages}")

    with st.expander("Lineage Distribution", expanded=False):
        lineage_counts = stats.get("lineage_distribution", {})
        if lineage_counts is not None and not lineage_counts.empty:
            lineage_counts.columns = ["Lineage", "Count"]
            st.dataframe(lineage_counts, hide_index=True, use_container_width=False)
        else:
            st.write("No lineage data available.")

def references(virus):
    st.subheader("References & Resources")
    refs_virus = REFERENCES.get(virus)
    refs_haplocov = REFERENCES.get(HAPLOCOV)
    refs = refs_virus + refs_haplocov
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

    haplocov()
    st.markdown("---")

    #haplocov_parameters(virus, config)
    #st.markdown("---")

    if virus in ["sars-cov-2"]: dataset_from_stats(virus, stats=df)
    else: dataset_from_df(virus, df=df)
    st.markdown("---")

    references(virus)
