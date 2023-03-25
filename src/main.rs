use bio::io::fasta;
use clap::{App, Arg};
use regex::Regex;
use rusqlite;


// read in in FASTA file and return a vector of sequences and metadata from header
fn parse_fasta(filename: &str) -> (Vec<(String, usize)>, Vec<(usize, usize, String, String, String, String, usize, usize)>) {
    let mut i: usize = 1; // protein number

    let mut seqs = Vec::new();
    let mut metadata = Vec::new();
    let reader = fasta::Reader::from_file(filename).unwrap();

    // regex to parse the header
    let re = Regex::new(r"(?x)
        (?:(sp|tr)\|)?
        \|(?P<protein_id>[^|]+)\|
        (\|[^|\s]+)?
        (?:\s(?P<protein_name>.+?))?\s
        (?:OS=(?P<species>.+?))?\s
        (?:OX=(?P<taxon_id>\d+))?\s
        (?:GN=(?P<gene>.+?))?\s
        (?:PE=(?P<pe_level>\d+))?\s
        (?:SV=(?P<sequence_version>\d+))?")
        .unwrap();

    for result in reader.records() {
        let record = result.unwrap();
        let seq_str = std::str::from_utf8(record.seq()).unwrap();
        seqs.push((seq_str.to_string(), i)); // store the sequence and the protein number
        
        // concatenate the id and description to get the full header
        let header = format!("{} {}", record.id(), record.desc().unwrap_or(""));

        // parse the header
        if let Some(caps) = re.captures(&header) {
            let protein_id = caps.name("protein_id").unwrap().as_str().to_string();
            let protein_name = caps.name("protein_name").map_or("".to_string(), |m| m.as_str().to_string());
            let species = caps.name("species").unwrap().as_str().to_string();
            let taxon_id: usize = caps.name("taxon_id").unwrap().as_str().parse().unwrap();
            let gene = caps.name("gene").map_or("".to_string(), |m| m.as_str().to_string());
            let pe_level: usize = caps.name("pe_level").unwrap().as_str().parse().unwrap();
            let sequence_version: usize = caps.name("sequence_version").unwrap().as_str().parse().unwrap();
            metadata.push((i, taxon_id, species, gene, protein_id, protein_name, pe_level, sequence_version));
        }
        i += 1;
    }

    (seqs, metadata)
}

// split the peptide into k-mers with a window size of 1 and store also the index of that k-mer
fn split_sequence(seq: &str, k: usize) -> Vec<(String, usize)> {
    let mut kmers = Vec::new();
    let mut i: usize = 0;
    while i + k <= seq.len() {
        kmers.push((seq[i..i + k].to_string(), i));
        i += 1;
    }
    kmers
}

// connect to SQLite DB, call it proteome.db
fn connect() -> rusqlite::Connection {
    rusqlite::Connection::open("proteome.db").unwrap()
}

// create a kmers --> index table in the DB
fn create_kmers_table(conn: &rusqlite::Connection) {
    conn.execute(
        "CREATE TABLE IF NOT EXISTS kmers (
            kmer             TEXT NOT NULL,
            idx              INTEGER NOT NULL
        )",
        rusqlite::params![],
    )
    .unwrap();
}

// create a protein metadata table in the DB
fn create_metadata_table(conn: &rusqlite::Connection) {
    conn.execute(
        "CREATE TABLE IF NOT EXISTS metadata (
            protein_number   INTEGER NOT NULL,
            taxon_id         INTEGER NOT NULL,
            species          TEXT NOT NULL,
            gene             TEXT NOT NULL,
            protein_id       TEXT NOT NULL,
            protein_name     TEXT NOT NULL,
            pe_level         INTEGER NOT NULL,
            sequence_version INTEGER NOT NULL
        )",
        rusqlite::params![],
    )
    .unwrap();
}

// insert kmers into the table
fn insert_kmers(conn: &mut rusqlite::Connection, kmers: &[(String, usize)], protein_count: &usize) {
    // Disable synchronous mode for faster bulk inserts
    conn.execute("PRAGMA synchronous = OFF", rusqlite::params![]).unwrap();

    let tx = conn.transaction().unwrap();
    let mut stmt = tx
        .prepare("INSERT INTO kmers (kmer, idx) VALUES (?1, ?2)")
        .unwrap();

    for kmer in kmers {
        stmt.execute(rusqlite::params![kmer.0, (protein_count * 1000000) + kmer.1])
            .unwrap();
    }

    drop(stmt); // Explicitly drop stmt before committing the transaction

    tx.commit().unwrap();

    // Re-enable synchronous mode
    conn.execute("PRAGMA synchronous = ON", rusqlite::params![]).unwrap();
}

// insert metadata into the table
fn insert_metadata(conn: &mut rusqlite::Connection, metadata: &[(usize, usize, String, String, String, String, usize, usize)]) {
    let tx = conn.transaction().unwrap();
    let mut stmt = tx
        .prepare("INSERT INTO metadata (protein_number, taxon_id, species, gene, protein_id, protein_name, pe_level, sequence_version) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)")
        .unwrap();

    for data in metadata {
        stmt.execute(rusqlite::params![data.0, data.1, data.2, data.3, data.4, data.5, data.6, data.7])
            .unwrap();
    }
    drop(stmt); // explicitly drop stmt before committing the transaction
    tx.commit().unwrap();
}

// create indices on the kmers and metadata tables
fn create_indices(conn: &mut rusqlite::Connection) {
    let tx = conn.transaction().unwrap();

    tx.execute("CREATE INDEX IF NOT EXISTS kmer_idx ON kmers (kmer)", rusqlite::params![])
        .unwrap();
    tx.execute("CREATE INDEX IF NOT EXISTS protein_number_idx ON metadata (protein_number)", rusqlite::params![])
        .unwrap();

    tx.commit().unwrap();
}

fn main() {
    let matches = App::new("Preprocess proteome.")
        .arg(
            Arg::with_name("proteome").short('p').long("proteome").value_name("FILE")
                .help("Input FASTA file").takes_value(true).required(true)
        )
        .arg(
            Arg::with_name("k").short('k').long("k_value").value_name("K")
                .help("Value of k for k-mers").takes_value(true).required(true),
        )
        .get_matches();

    let filename = matches.value_of("proteome").unwrap();
    let k: usize = matches.value_of("k").unwrap().parse()
        .unwrap_or_else(|_| {
            eprintln!("Error: k must be an integer");
            std::process::exit(1);
        });
    
    // parse proteome file and connect to DB
    let (seqs, metadata) = parse_fasta(filename);
    let mut conn = connect();

    // create metadata table and insert metadata
    create_metadata_table(&conn);
    insert_metadata(&mut conn, &metadata);

    // create kmers table and insert kmers
    create_kmers_table(&conn);
    for seq in seqs {
        let kmers = split_sequence(&seq.0, k);
        insert_kmers(&mut conn, &kmers, &seq.1);
    }

    // create indices
    create_indices(&mut conn);
}