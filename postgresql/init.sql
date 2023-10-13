use DATABASE db;

CREATE TABLE IF NOT EXISTS workspaces
(
    id INT NOT NULL,
    name VARCHAR(100) NOT NULL,
    data_dir VARCHAR(200) NOT NULL,
    created DATE NOT NULL,
    notes VARCHAR(10000) DEFAULT NULL,
    PRIMARY KEY (id)
);
