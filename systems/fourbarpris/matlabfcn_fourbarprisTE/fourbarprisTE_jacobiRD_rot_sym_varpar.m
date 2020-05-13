% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% JRD_rot [9x1]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = fourbarprisTE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_jacobiRD_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_jacobiRD_rot_sym_varpar: qJD has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbarprisTE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_jacobiRD_rot_sym_varpar: pkin has to be [3x1] (double)');
JRD_rot=NaN(9,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:49
	% EndTime: 2020-05-07 09:01:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:50
	% EndTime: 2020-05-07 09:01:50
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (568->25), mult. (226->29), div. (36->5), fcn. (20->2), ass. (0->21)
	t170 = qJ(1) + pkin(3);
	t178 = -pkin(2) + t170;
	t165 = pkin(1) + t178;
	t179 = -pkin(2) - t170;
	t166 = pkin(1) + t179;
	t184 = t165 - t166;
	t164 = pkin(1) - t178;
	t181 = t165 * t166;
	t183 = t184 * t164 + t181;
	t167 = 0.1e1 / t170;
	t168 = 0.1e1 / t170 ^ 2;
	t180 = qJD(1) / pkin(1);
	t177 = t164 * t181;
	t163 = pkin(1) - t179;
	t176 = t163 * t177;
	t160 = t183 * t163 - t177;
	t172 = sqrt(-t176);
	t169 = (t167 * t168);
	t161 = (-t167 + ((-pkin(1) ^ 2 + pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3)) * t169) + 0.2e1 * (t170 / 0.2e1 + pkin(3) / 0.2e1 + qJ(1) / 0.2e1) * t168) * t180;
	t158 = (t169 * t172 + (t167 * ((t164 - t184) * t163 + t183) / 0.2e1 + (-t168 / 0.2e1 + t167 * t160 / t176 / 0.8e1) * t160) / t172) * t180;
	t1 = [t161; -t158; 0; t158; t161; 0; 0; 0; 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:50
	% EndTime: 2020-05-07 09:01:50
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (568->25), mult. (226->29), div. (36->5), fcn. (20->2), ass. (0->21)
	t102 = qJ(1) + pkin(3);
	t110 = -pkin(2) + t102;
	t97 = pkin(1) + t110;
	t111 = -pkin(2) - t102;
	t98 = pkin(1) + t111;
	t116 = t97 - t98;
	t113 = t97 * t98;
	t96 = pkin(1) - t110;
	t115 = t116 * t96 + t113;
	t99 = 0.1e1 / t102;
	t100 = 0.1e1 / t102 ^ 2;
	t112 = qJD(1) / pkin(1);
	t109 = t96 * t113;
	t95 = pkin(1) - t111;
	t108 = t95 * t109;
	t92 = t115 * t95 - t109;
	t104 = sqrt(-t108);
	t101 = (t99 * t100);
	t93 = (t99 - ((-pkin(1) ^ 2 + pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3)) * t101) + 0.2e1 * (-t102 / 0.2e1 - qJ(1) / 0.2e1 - pkin(3) / 0.2e1) * t100) * t112;
	t90 = (-t101 * t104 + (-t99 * ((t96 - t116) * t95 + t115) / 0.2e1 + (t100 / 0.2e1 - t99 * t92 / t108 / 0.8e1) * t92) / t104) * t112;
	t1 = [t90; t93; 0; 0; 0; 0; -t93; t90; 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:01:49
	% EndTime: 2020-05-07 09:01:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (294->17), mult. (102->13), div. (12->2), fcn. (8->2), ass. (0->16)
	t133 = -pkin(3) - qJ(1);
	t130 = -pkin(2) - t133;
	t120 = pkin(1) + t130;
	t131 = -pkin(2) + t133;
	t121 = pkin(1) + t131;
	t136 = t120 - t121;
	t119 = pkin(1) - t130;
	t132 = t120 * t121;
	t135 = t136 * t119 + t132;
	t129 = t119 * t132;
	t128 = qJD(1) / pkin(1) / pkin(2);
	t118 = pkin(1) - t131;
	t127 = t118 * t129;
	t126 = t135 * t118 - t129;
	t117 = (-t126 ^ 2 / t127 / 0.8e1 - (t119 - t136) * t118 / 0.2e1 - t135 / 0.2e1) * (-t127) ^ (-0.1e1 / 0.2e1) * t128;
	t1 = [-t128; t117; 0; -t117; -t128; 0; 0; 0; 0;];
	JRD_rot = t1;
end