% Jacobian of implicit kinematic constraints of
% fourbar2IC
% with respect to passive joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% Phi_p [(no of constraints)x(no of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:37
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_p = fourbar2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:37:49
% EndTime: 2020-04-24 20:37:49
% DurationCPUTime: 0.04s
% Computational Cost: add. (4->3), mult. (4->4), div. (0->0), fcn. (4->4), ass. (0->2)
t7 = qJ(1) + qJ(2);
t1 = [-pkin(1) * sin(t7), -sin(qJ(3)) * pkin(2), 0; pkin(1) * cos(t7), cos(qJ(3)) * pkin(2), 0; -1, 1, -1;];
Phi_p = t1;
