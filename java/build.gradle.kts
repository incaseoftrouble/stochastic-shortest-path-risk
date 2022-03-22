plugins {
    java
    application
    idea
}

group = "ssp"
version = "0.2.0"

base {
    archivesName.set("ssp")
}

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17

    withSourcesJar()
    withJavadocJar()
}

dependencies {
    implementation(files("lib/gurobi.jar"))
    implementation("com.google.guava", "guava", "31.0.1-jre")
    implementation("com.google.code.gson", "gson", "2.8.9")
    implementation("info.picocli", "picocli", "4.6.3")
}

repositories {
    mavenCentral()
}

idea {
    module {
        isDownloadJavadoc = true
        isDownloadSources = true
    }
}

application {
    mainClass.set("ssp_risk.Main")
    applicationDefaultJvmArgs = listOf("-Xmx14G", "-Djava.util.logging.config.file=logging.properties")
}

tasks.test {
    useJUnitPlatform()
}