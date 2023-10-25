CREATE TABLE IF NOT EXISTS workspaces
(
    id INT NOT NULL,
    name VARCHAR(100) NOT NULL,
    data_dir VARCHAR(200) NOT NULL,
    created DATE NOT NULL,
    description VARCHAR(10000) DEFAULT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE IF NOT EXISTS adata
(
    work_id INT NOT NULL,
    id INT NOT NULL,
    name VARCHAR(100),
    adata BYTEA NOT NULL,
    notes VARCHAR(10000),
    created DATE NOT NULL,
    PRIMARY KEY (work_id, id),
    constraint fk_work_id foreign key (work_id) references workspaces (id)
)
