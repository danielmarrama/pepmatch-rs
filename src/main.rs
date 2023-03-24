use rusqlite;

// split the peptide into k-mers with a window size of 1 and store also the index of that k-mer
fn split_sequence(seq: &str, k: usize) -> Vec<(String, usize)> {
    let mut kmers = Vec::new();
    let mut i = 0;
    while i < seq.len() - k + 1 {
        kmers.push((seq[i..i + k].to_string(), i));
        i += 1;
    }
    kmers
}

// connect to sqlite db, call it proteome.db
fn connect() -> rusqlite::Connection {
    let conn = rusqlite::Connection::open("proteome.db").unwrap();
    conn
}

// create a kmers --> index table in the db
fn create_table(conn: &rusqlite::Connection) {
    conn.execute(
        "CREATE TABLE IF NOT EXISTS kmers (
                  kmer            TEXT NOT NULL,
                  idx             INTEGER NOT NULL
                  )",
        rusqlite::params![],
    )
    .unwrap();
}

// insert kmers into the table
fn insert_kmers(conn: &rusqlite::Connection, kmers: Vec<(String, usize)>) {
    let mut stmt = conn
        .prepare("INSERT INTO kmers (kmer, idx) VALUES (?1, ?2)")
        .unwrap();
    for kmer in kmers {
        stmt.execute(rusqlite::params![kmer.0, kmer.1])
            .unwrap();
    }
}

fn main() {
    let protein = "MLPGLALLLLAAWTARALEVPTDGNAGLLAEPQIAMFCGRLNMHMNVQNGKWDSDPSGTK";
    let kmers = split_sequence(protein, 3);
    println!("{:?}", kmers);
    let conn = connect();
    create_table(&conn);
    insert_kmers(&conn, kmers);
}

