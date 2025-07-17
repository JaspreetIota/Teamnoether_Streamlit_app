import streamlit as st
from Bio import Entrez
import pandas as pd
import re
from time import sleep
from io import BytesIO
import logging
import feedparser
import urllib.parse

# ---------- CONFIG ----------
Entrez.email = "your_email@example.com"  # Replace with your actual email
logging.basicConfig(level=logging.INFO)

# ---------- CONSTANTS ----------
PERSONAL_EMAIL_DOMAINS = ['gmail.com', 'yahoo.com', 'hotmail.com', 'outlook.com', 'aol.com']
BATCH_SIZE = 100

# ---------- FUNCTIONS ----------

def extract_email(text):
    match = re.search(r"[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+", text)
    return match.group(0).rstrip('.') if match else ""

def extract_university_name(affiliation_text):
    keywords = [
        'University', 'Institute', 'School', 'College', 'Hospital',
        'Laboratory', 'Lab', 'Centre', 'Center', 'Health', 'Academy'
    ]
    cleaned_aff = re.sub(r"[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+", "", affiliation_text)
    parts = [p.strip() for p in re.split(r'[;,]', cleaned_aff) if p.strip()]
    for part in reversed(parts):
        if any(k.lower() in part.lower() for k in keywords):
            if not re.search(r'\d', part) and '@' not in part and len(part.split()) > 1:
                return part
    return ""

def format_mla(authors, title, journal, volume, issue, year, pages, doi):
    author_str = ", ".join(authors)
    return f"{author_str}. \"{title}.\" *{journal}*, vol. {volume}, no. {issue}, {year}, pp. {pages}. doi:{doi}"

def get_google_news(query, max_articles=5):
    encoded_query = urllib.parse.quote(query)
    url = f'https://news.google.com/rss/search?q={encoded_query}'
    feed = feedparser.parse(url)
    news = []

    if not feed.entries:
        return []

    for entry in feed.entries[:max_articles]:
        news.append({
            'title': entry.title,
            'link': entry.link,
            'published': entry.published
        })
    return news

# ---------- STREAMLIT UI ----------
st.set_page_config(page_title="IOTA Tools", layout="wide")

menu = st.sidebar.selectbox("🔍 Select Tool", ["PubMed Article Extractor", "Google News Search"])

# ==========================
# 🚀 PubMed Article Extractor
# ==========================
if menu == "PubMed Article Extractor":
    st.title("🔬 IOTA's PubMed Article Extractor")

    search_term = st.text_input(
        "Enter your PubMed search term",
        '(Human Biology) AND ("united states"[Affiliation] OR USA[Affiliation]) AND (2022[Date - Publication])'
    )

    selected_countries = st.multiselect(
        "🌍 Select Countries (match in affiliation text)",
        options=[
            "USA", "United States", "United Kingdom", "Germany", "India", "Canada", "Australia",
            "France", "China", "Japan", "Brazil", "Italy", "Spain", "Netherlands", "Switzerland"
        ],
        default=["USA", "United States"]
    )

    retstart = st.number_input(
        "Start from record number (0 = first)",
        min_value=0, value=0, step=100
    )

    retmax = st.number_input(
        "How many records to fetch",
        min_value=10, max_value=10000, value=100, step=10
    )

    start_button = st.button("Fetch Articles")

    if not selected_countries:
        st.warning("⚠️ Please select at least one country to proceed.")

    if start_button and selected_countries:
        st.info(f"🔎 Searching PubMed (records {retstart} to {retstart + retmax - 1})...")
        try:
            search_handle = Entrez.esearch(
                db="pubmed",
                term=search_term,
                retstart=retstart,
                retmax=retmax
            )
            search_results = Entrez.read(search_handle)
        except Exception as e:
            st.error(f"❌ Failed to search PubMed: {e}")
            st.stop()

        pmids = search_results.get("IdList", [])
        if not pmids:
            st.error("No articles found for the given query and range.")
            st.stop()

        data = []
        progress = st.progress(0)

        for i, start in enumerate(range(0, len(pmids), BATCH_SIZE)):
            end = min(start + BATCH_SIZE, len(pmids))
            batch_pmids = pmids[start:end]

            try:
                fetch_handle = Entrez.efetch(db="pubmed", id=batch_pmids, rettype="xml")
                articles = Entrez.read(fetch_handle)
            except Exception as e:
                st.error(f"❌ Failed to fetch batch {start}-{end}: {e}")
                logging.exception("Fetch error")
                continue

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

                    elocations = article_data.get("ELocationID", [])
                    if isinstance(elocations, dict):
                        elocations = [elocations]
                    doi = "N/A"
                    for eloc in elocations:
                        if eloc.attributes.get("EIdType") == "doi":
                            doi = eloc.get("_", "N/A")
                            break

                    for author in article_data.get("AuthorList", []):
                        if "ForeName" in author and "LastName" in author:
                            full_name = f"{author['LastName']}, {author['ForeName']}"

                            for aff in author.get("AffiliationInfo", []):
                                aff_text = aff.get("Affiliation", "")
                                email = extract_email(aff_text)
                                if not email or any(email.lower().endswith(f"@{domain}") for domain in PERSONAL_EMAIL_DOMAINS):
                                    continue

                                matched_country = next(
                                    (country for country in selected_countries
                                     if re.search(rf'\b{re.escape(country)}\b', aff_text, re.IGNORECASE)),
                                    None
                                )
                                if not matched_country:
                                    continue

                                university = extract_university_name(aff_text)
                                if not university:
                                    continue

                                mla = format_mla([full_name], title, journal, volume, issue, year, pages, doi)

                                data.append({
                                    "PMID": pmid,
                                    "Author": full_name,
                                    "Email": email,
                                    "Country": matched_country,
                                    "University": university,
                                    "Affiliation": aff_text,
                                    "MLA Citation": mla
                                })
                                break
                except Exception as e:
                    logging.exception("Article processing error")
                    st.warning(f"⚠️ Skipped an article due to error: {e}")
                    continue

            sleep(0.5)
            progress.progress((i + 1) / ((len(pmids) - 1) // BATCH_SIZE + 1))

        if data:
            df = pd.DataFrame(data)
            st.success(f"✅ Completed! {len(data)} valid entries extracted.")
            st.dataframe(df)

            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name="Results")
            output.seek(0)

            st.download_button(
                label="📁 Download Excel",
                data=output,
                file_name="pubmed_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        else:
            st.error("❌ No valid data extracted. Try refining your search.")

# =====================
# 📰 Google News Search
# =====================
elif menu == "Google News Search":
    st.title("📰 Google News Article Finder")

    query = st.text_input("Enter a topic, company, or keyword", "BYD Auto")
    max_articles = st.slider("Number of articles to display", 1, 20, 5)

    if st.button("Search News"):
        with st.spinner("Fetching news..."):
            news_results = get_google_news(query, max_articles=max_articles)

        if news_results:
            st.success(f"✅ Found {len(news_results)} articles.")
            for i, article in enumerate(news_results, 1):
                st.markdown(f"**{i}. [{article['title']}]({article['link']})**")
                st.markdown(f"*Published:* {article['published']}\n")
        else:
            st.warning("⚠️ No news articles found or feed could not be loaded.")
