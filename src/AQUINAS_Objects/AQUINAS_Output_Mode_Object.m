classdef AQUINAS_Output_Mode_Object   < handle
    % AQUINAS_Output_Mode_Object An optional object used for defining the mode of output for AQUINAS
    % This Object is responsible for generating appropriate text files to present the results of the submitted analysis,
    % as well as determining whether AQUINAS GUI is to be called when the analysis ends.
    % AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.

    properties
        class_type = 'AQUINAS_Output_Mode_Object';
        mode                % The output mode that the User wants AQUINAS to run in.
                            % Appropriate values for this property are 'Silent', 'WriteFile', 'VarArgOut'
        level               % The detail level that the User wants information to be returned as output. Can be 'basic' or 'extensive'. Relevant for the 'VarArgOut'.
        file_path           % This has to be the absolute file path, including the desired file name but not the extension
        file_extension      % The type of the file where the output data will be saved.
                            % Appropriate values for this property are '.txt','.xls' and a '.<Custom>' for a User defined extension
        data                % Supported output data that can be selected for output in the desired file type.
                            % These have to be relevant to the type of analysis that is executed (see AQUINAS_Analysis_Object).
                            % The data that are available for each analysis type are the following:
                            %
        io_stream           % The output stream of this Output_Mode_Object
    end

    methods
        % Initiator method for an Output Mode Object
        function obj = AQUINAS_Output_Mode_Object(options)
            arguments
                options.mode { mustBeNonzeroLengthText } = 'Silent'
                options.level { mustBeNonzeroLengthText } = 'extensive'
                options.filePath { mustBeNonzeroLengthText } = ' '
                options.fileExtension { mustBeNonzeroLengthText } = '.aqn'
                options.data = []
            end
            obj.mode=options.mode; obj.file_path=options.filePath; obj.file_extension=options.fileExtension;
            obj.data=options.data; obj.io_stream=0;
            if strcmp(options.mode,'Silent') || strcmp(options.mode,'WriteFile') || strcmp(options.mode,'VarArgOut')
                obj.mode = options.mode;
            else
                error('AQUINAS Error : The only AQUINAS_Output_Mode_Object mode options supported are : Silent / WriteFile / VarArgOut.');
            end
            if strcmp(options.level,'basic') || strcmp(options.level,'extensive')
                obj.level = options.level;
            else
                error('AQUINAS Error : The only AQUINAS_Output_Mode_Object level options supported are either basic or extensive.');
            end
        end

        function obj = Open(obj)
           obj.io_stream=fopen(strcat(obj.file_extension,obj.file_path),'w');
        end

        function obj = Append(obj,analysis_data,analysis_obj)
            switch (analysis_obj.type)
                case ('LA')
                    for la_data=1:length(obj.data)

                    end
                case ('LBA')
                    line_output='';
                    for lba_data=1:size(obj.data,1)
                        switch obj.data(lba_data,:)
                            case ('LPF_crit')
                                line_output=strcat(line_output,num2str(analysis_data{1}{2}{1}),'    ');
                            case ('har_crit')
                        end
                    end
                    fprintf(obj.io_stream,'%s',strtrim(line_output));
            end
        end

        function obj =  Close(obj)
            fclose(obj.io_stream);
        end

        function obj = OpenWriteClose(obj,analysis_data,analysis_obj)
            obj.Open();
            obj.Append(analysis_data,analysis_obj);
            obj.Close();
        end

    end

    methods(Static)

        function CompleteOutput_PerSegmentHorizontally(io_fid, ANALYSIS, SEGMENTS, DataCell)

            charsPerVal = 17;

            fprintf(io_fid,'%s\n',      "oooooooooooooooooooooooooooo");
            fprintf(io_fid,'%s\n',      "<<<<<<<    AQUINAS   >>>>>>>");
            fprintf(io_fid,'%s\n\n',    "oooooooooooooooooooooooooooo");
            dt = fix(clock);
            fprintf(io_fid,'Date: %d/%d/%d   Time: %d:%d:%d\n\n\n',dt(3),dt(2),dt(1),dt(4),dt(5),dt(6));

            switch ANALYSIS.type
            case ('LA')
                % Main Headers
                fprintf(io_fid,'%s', " SEGMENT "); % 9 characters
                fprintf(io_fid,'%s', pad(pad("NODE",34,'left'),64,'right')); % 64 characters
                fprintf(io_fid,'%s',pad(pad("DEGREES OF FREEDOM",3*charsPerVal+9,'right'),6*charsPerVal,'left'));
                fprintf(io_fid,'%s',pad(pad("MEMBRANE STRESS RESULTANTS",ceil(1.5*charsPerVal)+13,'right'),3*charsPerVal,'left'));
                fprintf(io_fid,'%s',pad(pad("BENDING STRESS RESULTANTS",ceil(1.5*charsPerVal)+12,'right'),3*charsPerVal,'left'));
                fprintf(io_fid,'%s',pad(pad("STRAINS (MID-SURFACE)",ceil(1.5*charsPerVal)+10,'right'),3*charsPerVal,'left'));
                fprintf(io_fid,'%s',pad(pad("CURVATURES",ceil(1.5*charsPerVal)+5,'right'),3*charsPerVal,'left'));
                fprintf(io_fid,'%s',pad(pad("MERIDIONAL STRESSES",ceil(1.5*charsPerVal)+10,'right'),3*charsPerVal,'left'));
                fprintf(io_fid,'%s',pad(pad("CIRCUMFERENTIAL STRESSES",ceil(1.5*charsPerVal)+10,'right'),3*charsPerVal,'left'));
                fprintf(io_fid,'%s',pad(pad("VON MISES STRESSES",ceil(1.5*charsPerVal)+9,'right'),3*charsPerVal,'left'));

                % Sub-Headers
                fprintf(io_fid,'\n');
                fprintf(io_fid,'%s',"         "); % 9 characters
                fprintf(io_fid,'%s'," (GLOBAL ENUM) "); % 15 characters
                fprintf(io_fid,'%s'," (SEGMENT ENUM) "); % 16 characters
                fprintf(io_fid,'%s'," r coordinate "); % 15 characters
                fprintf(io_fid,'%s'," z coordinate "); % 15 characters
                fprintf(io_fid,'%s',pad(pad("u DOF",ceil(0.7*charsPerVal+3),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad("du/ds DOF",ceil(0.7*charsPerVal+4),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad("v DOF",ceil(0.7*charsPerVal+3),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad("dv/ds DOF",ceil(0.7*charsPerVal+4),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad("w DOF",ceil(0.7*charsPerVal+3),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad("dw/ds DOF",ceil(0.7*charsPerVal+4),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat('N',char(966)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat('N',char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat('N',char(966),char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat('M',char(966)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat('M',char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat('M',char(966),char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(949),char(966)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(949),char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(949),char(966),char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(954),char(966)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(954),char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(954),char(966),char(952)),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),char(966),',i'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),char(966),',m'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),char(966),',o'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),char(952),',i'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),char(952),',m'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),char(952),',o'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),'vM,i'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),'vM,m'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'%s',pad(pad(strcat(char(963),'vM,o'),ceil(0.7*charsPerVal+1),'left'),charsPerVal,'right'));
                fprintf(io_fid,'\n');
                node = 0;
                for S=1:length(SEGMENTS)
                    for Nd=1:(SEGMENTS(S).els+1)
                        fprintf(io_fid,'%s',pad(pad(num2str(S),5,'left'),9,'right')); % 9 characters
                        node = node + 1;
                        fprintf(io_fid,'%s',pad(pad(num2str(node),8,'left'),15,'right')); % 15 characters
                        fprintf(io_fid,'%s',pad(pad(num2str(Nd),9,'left'),16,'right')); % 16 characters
                        fprintf(io_fid,'%s',pad(num2str(SEGMENTS(S).r(Nd)),15,'left')); % 15 characters
                        fprintf(io_fid,'%s',pad(num2str(SEGMENTS(S).z(Nd)),15,'left')); % 15 characters
                        for i=1:6
                            entry = sprintf('%15.6e',DataCell{1}(6*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        for i=1:3
                            entry = sprintf('%15.6e',DataCell{4}(6*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        for i=4:6
                            entry = sprintf('%15.6e',DataCell{4}(6*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        for i=1:3
                            entry = sprintf('%15.6e',DataCell{2}(6*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        for i=4:6
                            entry = sprintf('%15.6e',DataCell{2}(6*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        for i=1:3:7
                            entry = sprintf('%15.6e',DataCell{5}(9*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        for i=2:3:8
                            entry = sprintf('%15.6e',DataCell{5}(9*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        for i=1:3
                            entry = sprintf('%15.6e',DataCell{6}(3*(node-1)+i));
                            fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); %charsPerVal characters
                        end
                        fprintf(io_fid,'\n');
                    end
                end

            case ('LBA')

                fprintf(io_fid,'%s\n\n',    "--- EIGENVALUES ---");
                fprintf(io_fid,'%s', " HARMONIC "); % 10 characters
                for N=1:ANALYSIS.no_eigenvalues
                    entry=sprintf(pad(strcat(" EIGENVALUE No ",num2str(N)," "),18,'left'));
                    fprintf(io_fid,'%s',entry); % 18 characters
                end
                fprintf(io_fid,'\n');
                for M=1:length(ANALYSIS.circumferential_modes)
                    fprintf(io_fid,pad(pad(num2str(ANALYSIS.circumferential_modes(M)),6,'left'),10,'right'));
                    for N=1:ANALYSIS.no_eigenvalues
                        entry = sprintf('%15.7f',DataCell{1}(ANALYSIS.no_eigenvalues*(M-1)+N));
                        fprintf(io_fid,'%s',pad(entry,18,'left'));
                    end
                    fprintf(io_fid,'\n');
                end

                fprintf(io_fid,'\n\n%s\n\n',    "--- EIGENMODES ---");
                fprintf(io_fid,'%s', " SEGMENT "); % 9 characters
                fprintf(io_fid,'%s', pad(pad("NODE",34,'left'),64,'right')); % 64 characters

                for M=1:length(ANALYSIS.circumferential_modes)
                    for N=1:ANALYSIS.no_eigenvalues
                        entry=sprintf('%s',strcat("DOF VALUES FOR CIRCUMFERENTIAL WAVENUMBER ",num2str(ANALYSIS.circumferential_modes(M))," AND EIGENMODE ",num2str(N))); % 28 characters + 3 characters for circumferential wavenumber value (circumferential wavenumber  up to 999) ceil(0.5*(28+3))=16
                        fprintf(io_fid,'%s', pad(pad(entry,0.5*4*charsPerVal+22,'left'),4*charsPerVal,'right')); % 4*charsPerVal characters
                    end
                end
                fprintf(io_fid,'\n');
                fprintf(io_fid,'%s',"         "); % 9 characters
                fprintf(io_fid,'%s'," (GLOBAL ENUM) "); % 15 characters
                fprintf(io_fid,'%s'," (SEGMENT ENUM) "); % 16 characters
                fprintf(io_fid,'%s'," r coordinate "); % 15 characters
                fprintf(io_fid,'%s'," z coordinate "); % 15 characters
                for M=1:length(ANALYSIS.circumferential_modes)
                    for N=1:ANALYSIS.no_eigenvalues
                        fprintf(io_fid,'%s',pad(pad("u DOF",ceil(0.5*charsPerVal),'left'),charsPerVal,'right')); % charsPerVal characters
                        fprintf(io_fid,'%s',pad(pad("v DOF",ceil(0.5*charsPerVal),'left'),charsPerVal,'right')); % charsPerVal characters
                        fprintf(io_fid,'%s',pad(pad("w DOF",ceil(0.5*charsPerVal),'left'),charsPerVal,'right')); % charsPerVal characters
                        fprintf(io_fid,'%s',pad(pad(strcat(char(946)," DOF"),ceil(0.5*charsPerVal),'left'),charsPerVal,'right')); % charsPerVal characters
                    end
                end
                fprintf(io_fid,'\n');

                node = 0;
                for S=1:length(SEGMENTS)
                    for Nd=1:(SEGMENTS(S).els+1)
                        fprintf(io_fid,'%s',pad(pad(num2str(S),5,'left'),9,'right')); % 9 characters
                        node = node + 1;
                        fprintf(io_fid,'%s',pad(pad(num2str(node),8,'left'),15,'right')); % 15 characters
                        fprintf(io_fid,'%s',pad(pad(num2str(Nd),9,'left'),16,'right')); % 16 characters
                        fprintf(io_fid,'%s',pad(num2str(SEGMENTS(S).r(Nd)),15,'left')); % 15 characters
                        fprintf(io_fid,'%s',pad(num2str(SEGMENTS(S).z(Nd)),15,'left')); % 15 characters
                        for M=1:length(ANALYSIS.circumferential_modes)
                            for N=1:ANALYSIS.no_eigenvalues
                                for i=1:4
                                    entry = sprintf('%15.6e',DataCell{1+ANALYSIS.no_eigenvalues*(M-1)+N}(4*(node-1)+i));
                                    fprintf(io_fid,'%s', pad(entry,charsPerVal,'left')); % charsPerVal characters
                                end
                            end
                        end
                        fprintf(io_fid,'\n');
                    end
                end
            end

        end

        function LA_struct = LA_structured_output(ANALYSIS, DOF_TRACKERS, SEGMENTS, DataCell, OStorage)

            LA_struct.Trackers = containers.Map('KeyType','char','ValueType','any');

            LA_struct.Shell_Geom.node   = cell([length(SEGMENTS) 1]);
            LA_struct.Shell_Geom.r = cell([length(SEGMENTS) 1]);
            LA_struct.Shell_Geom.z = cell([length(SEGMENTS) 1]);
            LA_struct.Shell_Geom.s = cell([length(SEGMENTS) 1]);
            LA_struct.DOFs.u = cell([length(SEGMENTS) 1]);
            LA_struct.DOFs.w = cell([length(SEGMENTS) 1]);
            LA_struct.DOFs.beta = cell([length(SEGMENTS) 1]);
            LA_struct.Stress_Reslts.Membrane.Nphi = cell([length(SEGMENTS) 1]);
            LA_struct.Stress_Reslts.Membrane.Ntheta = cell([length(SEGMENTS) 1]);
            LA_struct.Stress_Reslts.Membrane.Nphitheta = cell([length(SEGMENTS) 1]);
            LA_struct.Stress_Reslts.Bending.Mphi = cell([length(SEGMENTS) 1]);
            LA_struct.Stress_Reslts.Bending.Mtheta = cell([length(SEGMENTS) 1]);
            LA_struct.Stress_Reslts.Bending.Mphitheta = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_phi.i = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_phi.m = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_phi.o = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_theta.i = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_theta.m = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_theta.o = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_phitheta.m = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_vonMises.i = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_vonMises.m = cell([length(SEGMENTS) 1]);
            LA_struct.Stresses.Sig_vonMises.o = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Epsilon.Eps_phi.i = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Epsilon.Eps_phi.m = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Epsilon.Eps_phi.o = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Epsilon.Eps_theta.i = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Epsilon.Eps_theta.m = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Epsilon.Eps_theta.o = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Epsilon.Eps_phitheta.m = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Kappa.Kappa_phi = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Kappa.Kappa_theta = cell([length(SEGMENTS) 1]);
            LA_struct.Strains.Kappa.Kappa_phitheta = cell([length(SEGMENTS) 1]);

            for D=DOF_TRACKERS.keys

                TrackerDictTemp.ID         = DOF_TRACKERS(D{1}).ID;
                TrackerDictTemp.active_end = DOF_TRACKERS(D{1}).active_end;
                TrackerDictTemp.form       = DOF_TRACKERS(D{1}).form;
                TrackerDictTemp.DOF_type   = DOF_TRACKERS(D{1}).DOF_type;
                TrackerDictTemp.segID      = DOF_TRACKERS(D{1}).segID;
                TrackerDictTemp.global_DOF = DOF_TRACKERS(D{1}).global_DOF;
                TrackerDictTemp.DOF        = DOF_TRACKERS(D{1}).DOF;
                LA_struct.Trackers(DOF_TRACKERS(D{1}).name) = TrackerDictTemp;

            end

            el_counter = 0;

            for S=1:length(SEGMENTS)

                nod_top = SEGMENTS(S).top_dofIDs(4)/4;  top_3DOF = (nod_top-1)*3;   top_4DOF = (nod_top-1)*4;   top_6DOF = (nod_top-1)*6;   top_9DOF = (nod_top-1)*9;
                nod_bot = SEGMENTS(S).bot_dofIDs(4)/4;  bot_3DOF = nod_bot*3;       bot_4DOF = nod_bot*4;       bot_6DOF = nod_bot*6;       bot_9DOF = nod_bot*9;

                LA_struct.Shell_Geom.node{S}   = (nod_top : nod_bot)';
                LA_struct.Shell_Geom.r{S} = SEGMENTS(S).r';
                LA_struct.Shell_Geom.z{S} = SEGMENTS(S).z';
                LA_struct.Shell_Geom.s{S} = zeros(SEGMENTS(S).els,1);
                for i = 1:SEGMENTS(S).els
                    el_counter = el_counter +1;
                    LA_struct.Shell_Geom.s{S}(i+1) = LA_struct.Shell_Geom.s{S}(i) + 2*OStorage(11,el_counter);
                end
                LA_struct.DOFs.u{S} = DataCell{8}(top_4DOF+1:4:bot_4DOF);
                LA_struct.DOFs.w{S} = DataCell{8}(top_4DOF+3:4:bot_4DOF);
                LA_struct.DOFs.beta{S} = DataCell{8}(top_4DOF+4:4:bot_4DOF);
                LA_struct.Stress_Reslts.Membrane.Nphi{S} = DataCell{4}(top_6DOF+1:6:bot_6DOF);
                LA_struct.Stress_Reslts.Membrane.Ntheta{S} = DataCell{4}(top_6DOF+2:6:bot_6DOF);
                LA_struct.Stress_Reslts.Membrane.Nphitheta{S} = DataCell{4}(top_6DOF+3:6:bot_6DOF);
                LA_struct.Stress_Reslts.Bending.Mphi{S} = DataCell{4}(top_6DOF+4:6:bot_6DOF);
                LA_struct.Stress_Reslts.Bending.Mtheta{S} = DataCell{4}(top_6DOF+5:6:bot_6DOF);
                LA_struct.Stress_Reslts.Bending.Mphitheta{S} = DataCell{4}(top_6DOF+6:6:bot_6DOF);
                LA_struct.Stresses.Sig_phi.i{S} = DataCell{5}(top_9DOF+1:9:bot_9DOF);
                LA_struct.Stresses.Sig_phi.m{S} = DataCell{5}(top_9DOF+4:9:bot_9DOF);
                LA_struct.Stresses.Sig_phi.o{S} = DataCell{5}(top_9DOF+7:9:bot_9DOF);
                LA_struct.Stresses.Sig_theta.i{S} = DataCell{5}(top_9DOF+2:9:bot_9DOF);
                LA_struct.Stresses.Sig_theta.m{S} = DataCell{5}(top_9DOF+5:9:bot_9DOF);
                LA_struct.Stresses.Sig_theta.o{S} = DataCell{5}(top_9DOF+8:9:bot_9DOF);
                LA_struct.Stresses.Sig_phitheta.m{S} = DataCell{5}(top_9DOF+6:9:bot_9DOF);
                LA_struct.Stresses.Sig_vonMises.i{S} = DataCell{6}(top_3DOF+1:3:bot_3DOF);
                LA_struct.Stresses.Sig_vonMises.m{S} = DataCell{6}(top_3DOF+2:3:bot_3DOF);
                LA_struct.Stresses.Sig_vonMises.o{S} = DataCell{6}(top_3DOF+3:3:bot_3DOF);
                LA_struct.Strains.Epsilon.Eps_phi.i{S} = DataCell{3}(top_9DOF+1:9:bot_9DOF);
                LA_struct.Strains.Epsilon.Eps_phi.m{S} = DataCell{3}(top_9DOF+4:9:bot_9DOF);
                LA_struct.Strains.Epsilon.Eps_phi.o{S} = DataCell{3}(top_9DOF+7:9:bot_9DOF);
                LA_struct.Strains.Epsilon.Eps_theta.i{S} = DataCell{3}(top_9DOF+2:9:bot_9DOF);
                LA_struct.Strains.Epsilon.Eps_theta.m{S} = DataCell{3}(top_9DOF+5:9:bot_9DOF);
                LA_struct.Strains.Epsilon.Eps_theta.o{S} = DataCell{3}(top_9DOF+8:9:bot_9DOF);
                LA_struct.Strains.Epsilon.Eps_phitheta.m{S} = DataCell{3}(top_9DOF+3:9:bot_9DOF);
                LA_struct.Strains.Kappa.Kappa_phi{S} = DataCell{2}(top_6DOF+4:6:bot_6DOF);
                LA_struct.Strains.Kappa.Kappa_theta{S} = DataCell{2}(top_6DOF+5:6:bot_6DOF);
                LA_struct.Strains.Kappa.Kappa_phitheta{S} = DataCell{2}(top_6DOF+6:6:bot_6DOF);

            end

        end

        function LBA_struct = LBA_structured_output(ANALYSIS, SEGMENTS, DataCell, OStorage)

            LBA_struct.Shell_Geom.node   = cell([length(SEGMENTS) 1]);
            LBA_struct.Shell_Geom.r = cell([length(SEGMENTS) 1]);
            LBA_struct.Shell_Geom.z = cell([length(SEGMENTS) 1]);
            LBA_struct.Shell_Geom.s = cell([length(SEGMENTS) 1]);

            LBA_struct.CircWave = containers.Map('KeyType','int32','ValueType','any');

            el_counter = 0;

            for S=1:length(SEGMENTS)

                nod_top = SEGMENTS(S).top_dofIDs(4)/4;
                nod_bot = SEGMENTS(S).bot_dofIDs(4)/4;

                LBA_struct.Shell_Geom.node{S}   = nod_top : nod_bot;
                LBA_struct.Shell_Geom.r{S} = SEGMENTS(S).r;
                LBA_struct.Shell_Geom.z{S} = SEGMENTS(S).z;
                LBA_struct.Shell_Geom.s{S} = zeros(SEGMENTS(S).els,1);
                for i = 1:SEGMENTS(S).els
                    el_counter = el_counter +1;
                    LBA_struct.Shell_Geom.s{S}(i+1) = LBA_struct.Shell_Geom.s{S}(i) + 2*OStorage(11,el_counter);
                end

            end

            for C=1:length(ANALYSIS.circumferential_modes)

                CircWaveStructTemp.Circumferential_Wave_No = ANALYSIS.circumferential_modes(C);

                for M=1:ANALYSIS.no_eigenvalues

                    CircWaveStructTemp.EigenValue{M} = DataCell{1}(M,C);
                    CircWaveStructTemp.EigenMode{M}.u = cell([length(SEGMENTS) 1]);
                    CircWaveStructTemp.EigenMode{M}.v = cell([length(SEGMENTS) 1]);
                    CircWaveStructTemp.EigenMode{M}.w = cell([length(SEGMENTS) 1]);
                    CircWaveStructTemp.EigenMode{M}.beta = cell([length(SEGMENTS) 1]);

                    for S=1:length(SEGMENTS)

                        nod_top = SEGMENTS(S).top_dofIDs(4)/4;  top_4DOF = (nod_top-1)*4;
                        nod_bot = SEGMENTS(S).bot_dofIDs(4)/4;  bot_4DOF = nod_bot*4;

                        CircWaveStructTemp.EigenMode{M}.u{S} = DataCell{2}(top_4DOF+1:4:bot_4DOF,M,C);
                        CircWaveStructTemp.EigenMode{M}.v{S} = DataCell{2}(top_4DOF+2:4:bot_4DOF,M,C);
                        CircWaveStructTemp.EigenMode{M}.w{S} = DataCell{2}(top_4DOF+3:4:bot_4DOF,M,C);
                        CircWaveStructTemp.EigenMode{M}.beta{S} = DataCell{2}(top_4DOF+4:4:bot_4DOF,M,C);

                    end

                end
                LBA_struct.CircWave(ANALYSIS.circumferential_modes(C)) = CircWaveStructTemp;
            end

        end

        function GMNAstruct = GMNA_structured_output(output_level, ANALYSIS, SOLVER, SEGMENTS, DOF_TRACKERS, GMNAstep, OStorage)

            GMNAstruct.Shell_Geom.node   = cell([length(SEGMENTS) 1]);
            GMNAstruct.Shell_Geom.r = cell([length(SEGMENTS) 1]);
            GMNAstruct.Shell_Geom.z = cell([length(SEGMENTS) 1]);
            GMNAstruct.Shell_Geom.s = cell([length(SEGMENTS) 1]);

            el_counter = 0;
            Gauss = AQUINAS_Solver_Object.Gauss_Legendre_nodes_weights(SOLVER.No_Gauss_Stations);
            for S=1:length(SEGMENTS)
                nod_top = SEGMENTS(S).top_dofIDs(4)/4;
                nod_bot = SEGMENTS(S).bot_dofIDs(4)/4;
                GMNAstruct.Shell_Geom.node{S}   = nod_top : nod_bot;
                GMNAstruct.Shell_Geom.r{S} = SEGMENTS(S).r;
                GMNAstruct.Shell_Geom.z{S} = SEGMENTS(S).z;
                GMNAstruct.Shell_Geom.s{S} = zeros(SEGMENTS(S).els,1);
                for E = 1:SEGMENTS(S).els
                    el_counter = el_counter +1;
                    GMNAstruct.Shell_Geom.s{S}(E+1) = GMNAstruct.Shell_Geom.s{S}(E) + 2*OStorage(11,el_counter);
                    for I = 1:SOLVER.No_Gauss_Stations
                        eta = Gauss.nodes(I); eta2 = eta*eta; eta3 = eta2*eta; L = OStorage(11,el_counter);
                        r1 = GMNAstruct.Shell_Geom.r{S}(E); r2 = GMNAstruct.Shell_Geom.s{S}(E+1);
                        z1 = GMNAstruct.Shell_Geom.z{S}(E); z2 = GMNAstruct.Shell_Geom.z{S}(E+1);
                        drds1 = OStorage(2, el_counter); drds2 = OStorage(1, el_counter);
                        dzds1 = OStorage(4, el_counter); dzds2 = OStorage(3, el_counter);
                        N01 = (2-3*eta+eta3)*0.25; N11 = L*(1-eta-eta2+eta3)*0.25; N02 = (2+3*eta-eta3)*0.25; N12 = L*(-1-eta+eta2+eta3)*0.25; % Eq. 2 in Rotter & Teng (Vol. 31)

                        GMNAstruct.Shell_Geom.seg{S}.el{E}.gp{I}.s = GMNAstruct.Shell_Geom.s{S}(E) + L + Gauss.nodes(I)*L;
                        GMNAstruct.Shell_Geom.seg{S}.el{E}.gp{I}.r = N01*r1 + N02*r2 + N11*drds1 + N12*drds2; % Eqs. 1a & 2 in Rotter & Teng (Vol. 31)
                        GMNAstruct.Shell_Geom.seg{S}.el{E}.gp{I}.z = N01*z1 + N02*z2 + N11*dzds1 + N12*dzds2; % Eqs. 1b & 2 in Rotter & Teng (Vol. 31)
                    end
                end
            end

            GMNAstruct.Step = cell(length(GMNAstep),1); lsum = 0;
            for I = 1:length(GMNAstep)
                if strcmp(SOLVER.nonlinear_solver,'ArcLength'); lsum = lsum + GMNAstep{I}.l; GMNAstruct.Step{I}.l = lsum; end
                GMNAstruct.Step{I}.LPF = GMNAstep{I}.LPF;
                GMNAstruct.Step{I}.BuckledIntoMode = GMNAstep{I}.BuckledIntoMode;
                % Storing the results from bifurcation checks per GNA step
                GMNAstruct.Step{I}.BifurcationMode = containers.Map('KeyType','int32','ValueType','any'); % This is equivalent to LBA_struct.CircWave
                for Ckey = GMNAstep{I}.Bifurcation.keys
                    BifurcationStructTemp.Circumferential_Wave_No = Ckey{1};
                    for M=1:ANALYSIS.no_eigenvalues
                        BifurcationStructTemp.EigenValue{M} = GMNAstep{I}.Bifurcation(Ckey{1}).EigenValue(M);
                        BifurcationStructTemp.EigenMode{M}.u = cell(length(SEGMENTS),1);
                        BifurcationStructTemp.EigenMode{M}.v = cell(length(SEGMENTS),1);
                        BifurcationStructTemp.EigenMode{M}.w = cell(length(SEGMENTS),1);
                        BifurcationStructTemp.EigenMode{M}.beta = cell(length(SEGMENTS),1);
                        for S=1:length(SEGMENTS)
                            nod_top = SEGMENTS(S).top_dofIDs(4)/4;  top_4DOF = (nod_top-1)*4;
                            nod_bot = SEGMENTS(S).bot_dofIDs(4)/4;  bot_4DOF = nod_bot*4;
                            BifurcationStructTemp.EigenMode{M}.u{S} = GMNAstep{I}.Bifurcation(Ckey{1}).EigenVector(top_4DOF+1:4:bot_4DOF,M);
                            BifurcationStructTemp.EigenMode{M}.v{S} = GMNAstep{I}.Bifurcation(Ckey{1}).EigenVector(top_4DOF+2:4:bot_4DOF,M);
                            BifurcationStructTemp.EigenMode{M}.w{S} = GMNAstep{I}.Bifurcation(Ckey{1}).EigenVector(top_4DOF+3:4:bot_4DOF,M);
                            BifurcationStructTemp.EigenMode{M}.beta{S} = GMNAstep{I}.Bifurcation(Ckey{1}).EigenVector(top_4DOF+4:4:bot_4DOF,M);
                        end
                    end
                    GMNAstruct.Step{I}.BifurcationMode(Ckey{1}) = BifurcationStructTemp;
                end

                % Storing the nodal DOF values per step of the GMNA, as well as the strains and stresses for each Simpson station, for each Gauss point, for every element of every segment
                el_counter = 0;
                for S = 1:length(SEGMENTS)
                    nod_top = SEGMENTS(S).top_dofIDs(4)/4;  top_4DOF = (nod_top-1)*4;
                    nod_bot = SEGMENTS(S).bot_dofIDs(4)/4;  bot_4DOF = nod_bot*4;
                    GMNAstruct.Step{I}.DOFs.u{S} = GMNAstep{I}.Delta(top_4DOF+1:4:bot_4DOF);
                    GMNAstruct.Step{I}.DOFs.w{S} = GMNAstep{I}.Delta(top_4DOF+3:4:bot_4DOF);
                    GMNAstruct.Step{I}.DOFs.beta{S} = GMNAstep{I}.Delta(top_4DOF+4:4:bot_4DOF);
                    if strcmp(output_level,'extensive')
                        for E = 1:SEGMENTS(S).els
                            el_counter = el_counter + 1;
                            for J = 1:SOLVER.No_Gauss_Stations
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Stress_Reslts.Membrane.Nphi           = GMNAstep{I}.SIGMAS(1,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Stress_Reslts.Membrane.Ntheta         = GMNAstep{I}.SIGMAS(2,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Stress_Reslts.Membrane.Nphitheta      = GMNAstep{I}.SIGMAS(3,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Stress_Reslts.Bending.Mphi            = GMNAstep{I}.SIGMAS(4,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Stress_Reslts.Bending.Mtheta          = GMNAstep{I}.SIGMAS(5,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Stress_Reslts.Bending.Mphitheta       = GMNAstep{I}.SIGMAS(6,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Strains.Epsilon.Eps_phi               = GMNAstep{I}.epsilons(1,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Strains.Epsilon.Eps_theta             = GMNAstep{I}.epsilons(2,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Strains.Epsilon.Eps_phitheta          = GMNAstep{I}.epsilons(3,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Strains.Kappa.Kappa_phi               = GMNAstep{I}.epsilons(4,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Strains.Kappa.Kappa_theta             = GMNAstep{I}.epsilons(5,el_counter,J);
                                GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.Strains.Kappa.Kappa_phitheta          = GMNAstep{I}.epsilons(6,el_counter,J);
                                for K = 1:SOLVER.No_Simpson_Stations
                                    GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.sp{K}.Alpha                       = GMNAstep{I}.alphas(el_counter,J,K);
                                    GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.sp{K}.Stresses.Sig_phi            = GMNAstep{I}.sigmas(1,el_counter,J,K);
                                    GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.sp{K}.Stresses.Sig_theta          = GMNAstep{I}.sigmas(2,el_counter,J,K);
                                    GMNAstruct.Step{I}.seg{S}.el{E}.gp{J}.sp{K}.Stresses.Sig_phitheta       = GMNAstep{I}.sigmas(3,el_counter,J,K);
                                end
                            end
                        end
                    end
                end

                % Storing the nodal DOF values per step of the GNA, both in the condensed system (4 DOFs per node) and the original one (12 DOFs per element)
                GMNAstruct.Step{I}.DOF_Trackers = containers.Map;
                for Dkey=DOF_TRACKERS.keys
                    GMNAstruct.Step{I}.DOF_Trackers(Dkey{1}) = DOF_TRACKERS(Dkey{1}).DOF(I+1);
                end
            end
            GMNAstruct.Step{end}.termination_cause = GMNAstep{end}.termination_cause;
            if isfield(SOLVER.cip,'Rpl'); GMNAstep{end}.Rpl = SOLVER.cip.Rpl; end

        end

        % Method to visualise the progress of the geometrically and/or materially nonlinear analysis
        function [] = nonlinear_visualiser(ax,option,incr,iter,DELTA,norm_zeta_prev,zeta,alphas,sigmas,SOLVER,SEGMENTS,DOF_TRACKERS,OSTORAGE,include_nonlinear_M)

            switch option

            case 'deformedShape' % Deformed shape plot of current iteration

                nnodes = 0; for S = 1:length(SEGMENTS); nnodes = nnodes + SEGMENTS(S).els + 1; end
                nodcount = 0; Rcoord = zeros(nnodes,1); Zcoord = zeros(nnodes,1);
                for S=1:length(SEGMENTS)
                    Rcoord(1+nodcount:nodcount+SEGMENTS(S).els+1) = SEGMENTS(S).r';
                    Zcoord(1+nodcount:nodcount+SEGMENTS(S).els+1) = SEGMENTS(S).z';
                    nodcount = nodcount + SEGMENTS(S).els + 1;
                end
                plot(ax,Rcoord+DELTA(1:4:4*nnodes),Zcoord+DELTA(3:4:4*nnodes)); title(ax,strcat('Deformed Shape at Step--',num2str(incr)),strcat("for iteration ",num2str(iter))); xlabel(ax,'r coordinate [L]'); ylabel(ax,'z coordinate [L]'); grid(ax,'on');
                drawnow limitrate nocallbacks;

            case 'convergenceDevelopment' % Convergence criterion development

                if iter > 1 && ~isempty(zeta)
                    xlabel(ax,strcat("Iteration ",num2str(iter))); title(ax,'Convergence Criterion Development',strcat("for LPF = ",num2str(SOLVER.LPF)));
                    plot(ax,[iter-1 iter],[norm_zeta_prev norm(zeta,Inf)],'ro-');
                    plot(ax,[iter-1 iter],[SOLVER.zeta_t SOLVER.zeta_t],'b-');
                    plot(ax,[iter-1 iter],[SOLVER.zeta_d SOLVER.zeta_d],'b--');
                    legend(ax,{'change of nodal displacement magnitudes infinite norm','limit for convergence','limit for divergence'},'Location','best');
                    drawnow limitrate nocallbacks;
                end

            case 'translationalDOFtrackers'

                Dof_Tra_entries = {};
                for Dkey = DOF_TRACKERS.keys
                    if strcmp(DOF_TRACKERS(Dkey{1}).DOF_type,'u')||strcmp(DOF_TRACKERS(Dkey{1}).DOF_type,'w')
                        Dof_Tra_entries(end+1) = Dkey;
                    end
                end
                for Dkey = DOF_TRACKERS.keys
                    if strcmp(DOF_TRACKERS(Dkey{1}).DOF_type,'u')||strcmp(DOF_TRACKERS(Dkey{1}).DOF_type,'w')
                        plot(ax,[DOF_TRACKERS(Dkey{1}).DOF(incr) DOF_TRACKERS(Dkey{1}).DOF(incr+1)],[DOF_TRACKERS(Dkey{1}).LPF(incr) DOF_TRACKERS(Dkey{1}).LPF(incr+1)],'Color',DOF_TRACKERS(Dkey{1}).colour,'Marker',DOF_TRACKERS(Dkey{1}).marker,'MarkerSize',DOF_TRACKERS(Dkey{1}).markerSize);
                    end
                end
                legend(ax,Dof_Tra_entries,'Location','best');
                drawnow limitrate nocallbacks;

            case 'rotationalDOFtrackers'

                Dof_Rot_entries = {};
                for Dkey = DOF_TRACKERS.keys
                    if strcmp(DOF_TRACKERS(Dkey{1}).DOF_type,'beta') || strcmp(DOF_TRACKERS(Dkey{1}).DOF_type,'b')
                        Dof_Rot_entries(end+1) = Dkey;
                    end
                end
                for Dkey = DOF_TRACKERS.keys
                    if strcmp(DOF_TRACKERS(Dkey{1}).DOF_type,'beta')
                        plot(ax,[DOF_TRACKERS(Dkey{1}).DOF(incr) DOF_TRACKERS(Dkey{1}).DOF(incr+1)],[DOF_TRACKERS(Dkey{1}).LPF(incr) DOF_TRACKERS(Dkey{1}).LPF(incr+1)],'Color',DOF_TRACKERS(Dkey{1}).colour,'Marker',DOF_TRACKERS(Dkey{1}).marker);
                    end
                end
                legend(ax,Dof_Rot_entries,'Location','best');

            case 'alphas'

                numelem = size(alphas,1);
                numGauss = size(alphas,2);
                numSimpson = size(alphas,3);
                alphas2D = zeros(numSimpson,numelem*numGauss);
                Scoord = zeros(numelem*numGauss,1);
                Gauss = AQUINAS_Solver_Object.Gauss_Legendre_nodes_weights(numGauss);
                Stotal = 0;

                for E = 1:numelem
                    Stotal = Stotal + 2*OSTORAGE(11,E);
                    for I = 1:numGauss
                        Scoord((E-1)*numGauss + I) = OSTORAGE(11,E) + Gauss.nodes(I)*OSTORAGE(11,E);
                        if E > 1; Scoord((E-1)*numGauss + I) = Scoord((E-1)*numGauss + I) + 2*sum(OSTORAGE(11,1:(E-1))); end
                        for J = 1:numSimpson
                            alphas2D(J,numGauss*(E-1)+I) = alphas(E,I,J);
                        end
                    end
                end
                Scoord = Scoord/Stotal;
                if any(any(alphas2D))
                    [X,Y] = meshgrid(Scoord,1:numSimpson);
                    contourf(ax,X,Y,alphas2D,1); cmap = [1.0 1.0 1.0; 1.0 0.4 0.0]; colormap(ax,cmap); colorbar(ax);
                    title(ax,strcat("Yielded points during iteration ",num2str(iter)," of load step ",num2str(incr)," for LPF = ",num2str(SOLVER.LPF))); xlabel(ax,'Normalised arc length coordinate along the meridian of the shell (all elements)'); ylabel(ax,'Through thickness integration station');
                end

            case 'sigmaVonMises'

                numelem = size(sigmas,2);
                numGauss = size(sigmas,3);
                numSimpson = size(sigmas,4);
                sigmabars2D = zeros(numSimpson,numelem*numGauss);
                Scoord = zeros(numelem*numGauss,1);
                Gauss = AQUINAS_Solver_Object.Gauss_Legendre_nodes_weights(numGauss);
                Stotal = 0;

                for E = 1:numelem
                    Stotal = Stotal + 2*OSTORAGE(11,E);
                    for I = 1:numGauss
                        Scoord((E-1)*numGauss + I) = OSTORAGE(11,E) + Gauss.nodes(I)*OSTORAGE(11,E);
                        if E > 1; Scoord((E-1)*numGauss + I) = Scoord((E-1)*numGauss + I) + 2*sum(OSTORAGE(11,1:(E-1))); end
                        for J = 1:numSimpson
                            sigmabars2D(J,numGauss*(E-1)+I) = sqrt(sigmas(1,E,I,J)^2 + sigmas(2,E,I,J)^2 - sigmas(1,E,I,J)*sigmas(2,E,I,J));
                        end
                    end
                end
                Scoord = Scoord/Stotal;
                [X,Y] = meshgrid(Scoord,1:numSimpson);
                contourf(ax,X,Y,sigmabars2D,20); colormap(ax,jet); colorbar(ax);
                title(ax,strcat("Von Mises stresses during iteration ",num2str(iter)," of load step ",num2str(incr)," for LPF = ",num2str(SOLVER.LPF))); xlabel(ax,'Normalised arc length coordinate along the meridian of the shell (all elements)'); ylabel(ax,'Through thickness integration station');

            case 'sigmaphis'

                numelem = size(sigmas,2);
                numGauss = size(sigmas,3);
                numSimpson = size(sigmas,4);
                sigmaphis2D = zeros(numSimpson,numelem*numGauss);
                Scoord = zeros(numelem*numGauss,1);
                Gauss = AQUINAS_Solver_Object.Gauss_Legendre_nodes_weights(numGauss);
                Stotal = 0;

                for E = 1:numelem
                    Stotal = Stotal + 2*OSTORAGE(11,E);
                    for I = 1:numGauss
                        Scoord((E-1)*numGauss + I) = OSTORAGE(11,E) + Gauss.nodes(I)*OSTORAGE(11,E);
                        if E > 1; Scoord((E-1)*numGauss + I) = Scoord((E-1)*numGauss + I) + 2*sum(OSTORAGE(11,1:(E-1))); end
                        for J = 1:numSimpson
                            sigmaphis2D(J,numGauss*(E-1)+I) = sigmas(1,E,I,J);
                        end
                    end
                end
                Scoord = Scoord/Stotal;
                [X,Y] = meshgrid(Scoord,1:numSimpson);
                contourf(ax,X,Y,sigmaphis2D,20); colormap(ax,jet); colorbar(ax);
                title(ax,strcat("Meridional stresses during iteration ",num2str(iter)," of load step ",num2str(incr)," for LPF = ",num2str(SOLVER.LPF))); xlabel(ax,'Normalised arc length coordinate along the meridian of the shell (all elements)'); ylabel(ax,'Through thickness integration station');

            case 'sigmathetas'

                numelem = size(sigmas,2);
                numGauss = size(sigmas,3);
                numSimpson = size(sigmas,4);
                sigmathetas2D = zeros(numSimpson,numelem*numGauss);
                Scoord = zeros(numelem*numGauss,1);
                Gauss = AQUINAS_Solver_Object.Gauss_Legendre_nodes_weights(numGauss);
                Stotal = 0;

                for E = 1:numelem
                    Stotal = Stotal + 2*OSTORAGE(11,E);
                    for I = 1:numGauss
                        Scoord((E-1)*numGauss + I) = OSTORAGE(11,E) + Gauss.nodes(I)*OSTORAGE(11,E);
                        if E > 1; Scoord((E-1)*numGauss + I) = Scoord((E-1)*numGauss + I) + 2*sum(OSTORAGE(11,1:(E-1))); end
                        for J = 1:numSimpson
                            sigmathetas2D(J,numGauss*(E-1)+I) = sigmas(2,E,I,J);
                        end
                    end
                end
                Scoord = Scoord/Stotal;
                [X,Y] = meshgrid(Scoord,1:numSimpson);
                contourf(ax,X,Y,sigmathetas2D,20); colormap(ax,jet); colorbar(ax);
                title(ax,strcat("Circumferential stresses during iteration ",num2str(iter)," of load step ",num2str(incr)," for LPF = ",num2str(SOLVER.LPF))); xlabel(ax,'Normalised arc length coordinate along the meridian of the shell (all elements)'); ylabel(ax,'Through thickness integration station');

            end

        end

    end

end


% BSD 3-Clause License
%  
% Copyright (c) 2023, Mr Achilleas Filippidis and Dr Adam Jan Sadowski of 
% Imperial College London. All rights reserved.
%  
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%  
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%  
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  
% 3. Neither the name of the copyright holder, nor of Imperial College, nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%  
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%  
% THE USE OF THE MATLAB LANGUAGE DOES NOT IMPLY ENDORSEMENT BY MATHWORKS.s