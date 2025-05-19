import streamlit as st
from Bio import Entrez
import pandas as pd
import re
from time import sleep
from io import BytesIO  # Added for Excel download

# ---------- CONFIG ----------
Entrez.email = "your_email@example.com"

def extract_university_name(affiliation_text):
    keywords = ['University', ' Institute', 'School', 'College', 'Hospital',
                'Laboratory', 'Lab', 'Centre', 'Center', 'Health', 'Academy']
    cleaned_aff = re.sub(r"[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+", "", affiliation_text)
    parts = [p.strip() for p in re.split(r'[;,]', cleaned_aff) if p.strip()]
    for part in reversed(parts):
        if any(k.lower() in part.lower() for k in keywords):
            if not re.search(r'\d', part) and '@' not in part and len(part.split()) > 1:
                return part
    return ""

def extract_email(text):
    match = re.search(r"[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+", text)
    return match.group(0).rstrip('.') if match else ""

def format_mla(authors, title, journal, volume, issue, year, pages, doi):
    author_str = ", ".join(authors)
    return f"{author_str}. \"{title}.\" *{journal}*, vol. {volume}, no. {issue}, {year}, pp. {pages}. doi:{doi}"

# ---------- STREAMLIT UI ----------
st.title("🔬 PubMed Article Extractor")
search_term = st.text_input("Enter your PubMed search term", '(Human Biology) AND ("united states"[Affiliation] OR USA[Affiliation]) AND (2022[Date - Publication])')
max_results = st.number_input("Max Results", min_value=10, max_value=10000, value=100)
start_button = st.button("Fetch Articles")

if start_button:
    st.info("Searching PubMed...")
    search_handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_results)
    search_results = Entrez.read(search_handle)
    pmids = search_results["IdList"]

    BATCH_SIZE = 100
    data = []

    for start in range(0, len(pmids), BATCH_SIZE):
        end = min(start + BATCH_SIZE, len(pmids))
        batch_pmids = pmids[start:end]

        fetch_handle = Entrez.efetch(db="pubmed", id=batch_pmids, rettype="xml")
        articles = Entrez.read(fetch_handle)

        for article in articles['PubmedArticle']:
            try:
                pmid = str(article['MedlineCitation']['PMID'])
                article_data = article['MedlineCitation']['Article']
                title = article_data.get("ArticleTitle", "No Title")
                journal = article_data["Journal"]["Title"]
                volume = article_data["Journal"]["JournalIssue"].get("Volume", "N/A")
                issue = article_data["Journal"]["JournalIssue"].get("Issue", "N/A")
                year = article_data["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
                pages = article_data.get("Pagination", {}).get("MedlinePgn", "N/A")

                doi = "N/A"
                for eloc in article_data.get("ELocationID", []):
                    if isinstance(eloc, dict) and eloc.attributes.get("EIdType") == "doi":
                        doi = eloc.get("_", "N/A")

                for author in article_data.get("AuthorList", []):
                    if "ForeName" in author and "LastName" in author:
                        full_name = f"{author['LastName']}, {author['ForeName']}"
                        for aff in author.get("AffiliationInfo", []):
                            aff_text = aff.get("Affiliation", "")
                            email = extract_email(aff_text)
                            if not email or "gmail.com" in email:
                                continue

                            country_match = re.search(r'\b(USA|United States)\b', aff_text, re.IGNORECASE)
                            if not country_match:
                                continue

                            university = extract_university_name(aff_text)
                            if not university:
                                continue

                            mla = format_mla([full_name], title, journal, volume, issue, year, pages, doi)

                            data.append({
                                "PMID": pmid,
                                "Author": full_name,
                                "Email": email,
                                "Country": country_match.group(0),
                                "University": university,
                                "Affiliation": aff_text,
                                "MLA Citation": mla
                            })
                            break
            except Exception as e:
                st.warning(f"Error: {e}")
                continue

        sleep(0.5)
        st.success(f"✅ Processed {end}/{len(pmids)} articles...")

    if data:
        df = pd.DataFrame(data)
        st.dataframe(df)

        # ✅ Fix: Create Excel file in memory
        output = BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False, sheet_name="Results")
        output.seek(0)

        # ✅ Working download button
        st.download_button(
            label="📁 Download Excel",
            data=output,
            file_name="pubmed_results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    else:
        st.error("No valid results found.")
