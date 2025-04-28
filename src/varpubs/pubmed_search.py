from pathlib import Path
from typing import List
from sqlmodel import Session
from Bio import Entrez
from pubmed_db import PubmedDB, PubmedArticle  
import tkinter as tk
from tkinter import ttk
import duckdb

while True:
    email = input("Please enter your email address for NCBI Entrez access: ").strip()
    if "@" in email and "." in email:
        Entrez.email = email
        break
    else:
        print("Invalid email address. Please try again.")

db = PubmedDB(path=Path("pubmed_articles.db"))
db.deploy()


def search_and_fetch_articles(mutations: List[str], max_results_per_mutation: int = 100) -> List[PubmedArticle]:
    all_articles = []

    for mutation in mutations:
        print(f"Searching for '{mutation}' in PubMed...")
        articles = []
        search_handle = Entrez.esearch(db="pubmed", term=mutation, retmax=max_results_per_mutation)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        id_list = search_results.get("IdList", [])
        if not id_list:
            print(f"No articles found for mutation: {mutation}")
            continue
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="abstract", retmode="xml")
        fetch_data = Entrez.read(fetch_handle)
        fetch_handle.close()
        for article in fetch_data['PubmedArticle']:
            try:
                pmid = int(article['MedlineCitation']['PMID'])
                article_title = article['MedlineCitation']['Article']['ArticleTitle']
                abstract_parts = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                abstract = " ".join(abstract_parts) if abstract_parts else "No abstract available"
                authors = []
                for author in article['MedlineCitation']['Article'].get('AuthorList', []):
                    if "LastName" in author and "ForeName" in author:
                        authors.append(f"{author['ForeName']} {author['LastName']}")
                authors_str = ", ".join(authors) if authors else "Unknown authors"
                journal = article['MedlineCitation']['Article']['Journal']['Title']
                pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', "Unknown year")
                article_ids = article.get('PubmedData', {}).get('ArticleIdList', [])
                doi = "No DOI"
                for id_item in article_ids:
                    if id_item.attributes.get('IdType') == 'doi':
                        doi = str(id_item)
                pubmed_article = PubmedArticle(
                    pmid=pmid,
                    title=article_title,
                    abstract=abstract,
                    authors=authors_str,
                    journal=journal,
                    pub_date=pub_date,
                    doi=doi
                )
                articles.append(pubmed_article)
            except Exception as e:
                print(f"Skipping article due to error: {e}")

        print(f"Found {len(articles)} articles for '{mutation}'")
        all_articles.extend(articles)
    return all_articles

def insert_articles_into_db(db: PubmedDB, articles: list[PubmedArticle]):
    from sqlmodel import Session

    if not articles:
        return  

    with Session(db.engine, expire_on_commit=False) as session:
        values_placeholders = ", ".join(["(?, ?, ?, ?, ?, ?, ?)"] * len(articles))
        values = []
        for article in articles:
            values.extend([
                article.pmid,
                article.title,
                article.abstract,
                article.authors,
                article.journal,
                article.pub_date,
                article.doi,
            ])
        sql = f"""
        INSERT OR REPLACE INTO pubmedarticle (pmid, title, abstract, authors, journal, pub_date, doi)
        VALUES {values_placeholders}
        """
        for article in articles:
            exists = session.query(PubmedArticle).filter_by(pmid=article.pmid).first()
            if exists:
                continue
            session.add(article)
        session.commit()


user_input = input("Enter mutations separated by commas (e.g. BRCA1, TP53, KRAS): ")

mutations = [mutation.strip() for mutation in user_input.split(",")]
print("Mutations to search for:", mutations)

found_articles = search_and_fetch_articles(mutations, max_results_per_mutation=100)

if found_articles:
    insert_articles_into_db(db, found_articles)
    print("Articles inserted:")
    for article in found_articles:
        print(f"- {article.title} (PMID: {article.pmid}) [{article.journal}]")
else:
    print(" No articles were found.")


def show_db_in_tkinter(db_path, mutations):
    conn = duckdb.connect(database=str(db_path), read_only=True)
    cursor = conn.cursor()

    root = tk.Tk()
    root.title("PubMed Articles by Mutation")

    notebook = ttk.Notebook(root)
    notebook.pack(expand=True, fill="both")

    for mutation in mutations:
        frame = ttk.Frame(notebook)
        notebook.add(frame, text=mutation)

        table_frame = tk.Frame(frame)  
        table_frame.pack(expand=True, fill="both")

        columns = ("pmid", "title", "authors", "journal", "pub_date", "doi", "abstract")
        tree = ttk.Treeview(table_frame, columns=columns, show="headings")

        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(table_frame, orient="horizontal", command=tree.xview)

        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        table_frame.grid_rowconfigure(0, weight=1)
        table_frame.grid_columnconfigure(0, weight=1)

        for col in columns:
            tree.heading(col, text=col)
            if col == "abstract":
                tree.column(col, width=120, stretch=False) 
            else:
                tree.column(col, width=400, stretch=True)   

        cursor.execute("""
            SELECT pmid, title, authors, journal, pub_date, doi, abstract
            FROM pubmedarticle
            WHERE title ILIKE ? OR abstract ILIKE ?
        """, (f"%{mutation}%", f"%{mutation}%"))

        rows = cursor.fetchall()

        for row in rows:
            short_row = list(row)
            short_row[-1] = "Show Abstract"
            tree.insert("", tk.END, values=short_row)

        def on_item_double_click(event, mutation=mutation):
            selected_item = tree.selection()
            if not selected_item:
                return
            item = tree.item(selected_item)
            values = item["values"]
            pmid = values[0]

            cursor.execute("SELECT abstract FROM pubmedarticle WHERE pmid = ?", (pmid,))
            result = cursor.fetchone()
            if result:
                abstract = result[0]
                popup = tk.Toplevel(root)
                popup.title(f"Abstract for PMID {pmid} ({mutation})")
                text = tk.Text(popup, wrap="word", width=100, height=30)
                text.insert("1.0", abstract)
                text.config(state="disabled")
                text.pack(expand=True, fill="both")

        tree.bind("<Double-1>", on_item_double_click)


    root.mainloop()



db_file = Path("pubmed_articles.db")
db.engine.dispose()
show_db_in_tkinter(db_file, mutations)