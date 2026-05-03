FROM nginxinc/nginx-unprivileged:alpine
COPY nginx.conf /etc/nginx/conf.d/default.conf
COPY pn-junction-simulator.html /usr/share/nginx/html/
EXPOSE 8080
