use std::io::{BufWriter, Write};
use chrono::Utc;

pub fn log_event(log_stream: &mut Option<&mut BufWriter<Box<dyn Write>>>, msg: &str) {
    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\t{msg}", now.to_rfc2822())
            .expect("E: Failed to write log");
    }
}